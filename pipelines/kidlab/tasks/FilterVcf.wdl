version 1.0

## Copyright Marco Baggio, KU Leuven, 2023
## For questions: kuleuven@mawumag.com
## 
## This WDL defines tasks used for annotation of WES/WGS data aligned with the GATK pipeline.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

task FilterVcf {
  input {
    File input_vcf
    File ped_file
    String code
  }


  command <<<
    input=~{input_vcf}
    code=~{code}
    sex=$(awk '{if ($1=="'"$code"'" && $6=="2") print $5}' ~{ped_file})
    num=$(awk '{if ($1=="'"$code"'") print}' ~{ped_file} | wc -l)
    if (( num==1 )); then
      ped="FAM proband 0 0 $sex 2"
    elif (( num==2 )); then
      ped="FAM mother 0 0 2 1\nFAM proband 0 mother $sex 2"
    elif (( num>=3 )); then
      ped="FAM father 0 0 1 1\nFAM mother 0 0 2 1\nFAM proband father mother $sex 2"
    else
      exit 1
    fi

    if [ ! -f "$input" ]; then
      echo "File $input not found"
      exit 1
    fi

    format='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%AS_VQSLOD[\t%GT][\t%DP][\t%AD]\t%GeneticModels\t'
    csq=$(bcftools +split-vep -l "$input" | cut -f 2 | tr '\n' '\t' | sed 's/\t$/\\n/' | sed 's/\t/\\t%/g' | sed 's/^/%/') 
    format+=$csq

    proband=0
    father=1
    mother=2

    samples=($(bcftools query -l GC115416.vcf.gz))

    for i in $(seq 0 $((${#samples[@]}-1))); do
      if [[ ${samples[$i]} == "proband" ]]; then
        proband=$i
      elif [[ ${samples[$i]} == "father" ]]; then
        father=$i
      elif [[ ${samples[$i]} == "mother" ]]; then
        mother=$i
      fi
    done

    # Preliminary splitting and filtration

    bcftools view -i'GT['"$proband"']="alt" && FMT/DP['"$proband"']>9' "$input" | \
    bcftools +split-vep -c - -dx | \
    sed 's/ID=gnomAD_nhomalt,Number=.,Type=String/ID=gnomAD_nhomalt,Number=1,Type=Integer/' | \
    sed "s/^chrX/X/" | \
    genmod models -f <(echo -e "$ped") -k "Gene" - | \
    sed "s/^X/chrX/" > temp.vcf

    bcftools query -H -f"$format" temp.vcf > "$code".anno.txt

    bcftools view -f PASS temp.vcf | \
    bcftools view -i'AF<0.1' | \
    bcftools +split-vep -x -s all:5_prime_utr+ | \
    bcftools view -e'(CLIN_SIG~"benign" & CLIN_SIG!~"pathogenic") | gnomAD_nhomalt > 5' > temp2.vcf

    bcftools view -e'MAX_AF>0.02 | INFO/PID="."' temp2.vcf | \
    bcftools query -H -f"$format" > "$code".PID.txt

    bcftools view -e'MAX_AF>0.02 | INFO/PIDg="."' temp2.vcf | \
    bcftools query -H -f"$format" > "$code".PID_extra.txt

    bcftools view -e'MAX_AF>0.02 | INFO/Connectome="." | INFO/PID="1" | INFO/PIDg="1" | (CADD_PHRED < 15 & MSC=".") | (CADD_PHRED < MSC) | (GT['"$proband"']="het" & Domino<0.5)' temp2.vcf | \
    bcftools query -H -f"$format" > "$code".connectome.txt

    bcftools view -e'MAX_AF>0.02 | INFO/Connectome!="." | INFO/PID="1" | INFO/PIDg="1" | (CADD_PHRED < 15 & MSC=".") | (CADD_PHRED < MSC) | (GT['"$proband"']="het" & Domino<0.5)' temp2.vcf | \
    bcftools query -H -f"$format" > "$code".rest.txt

    bcftools view -e'MAX_AF>0.02 | AD['"$proband"':1]<0.25 * FORMAT/DP['"$proband"'] | AD['"$proband"':1]>0.75 * FORMAT/DP['"$proband"']' temp2.vcf | \
    bcftools view -i'Existing_variation="." & gnomAD="."' | \
    bcftools query -H -f"$format" > "$code".private.txt

    if (( num >= 3 )); then
      bcftools view -i'GT['"$father"']="ref" && GT['"$mother"']="ref"' temp2.vcf | \
      bcftools query -H -f"$format" > "$code".AD_denovo.txt

      bcftools view -i'(GeneticModels~"AR_comp")' temp2.vcf | \
      bcftools query -H -f"$format" > "$code".AR_comp.txt

      bcftools view -e'MAX_AF>0.141' temp2.vcf | \
      bcftools view -i'GT['"$proband"']="hom" && GT['"$father"']="het" && GT['"$mother"']="het"' | \
      bcftools query -H -f"$format" > "$code".AR_hom.txt
    fi

    for filter in "$code".*; do
      sed -i -r '1 s/\[[0-9]+\]//g' "$filter"
      sed -i -r '1 s/^#[ ]?//' "$filter"
    done

    zip "$code".anno.txt.zip "$code".anno.txt && rm "$code".anno.txt

    python3 << EOF
    import pandas as pd
    import csv
    from xlsxwriter.workbook import Workbook
    import os

    if $num>=3:
      df=pd.read_table("$code.AR_comp.txt",dtype=str)
      genes=df['Gene'].unique()
      for gene in genes:
        max_af = pd.to_numeric(df.loc[df['Gene']==gene]['MAX_AF'],errors='coerce').replace('NaN',0).max()
        min_af = pd.to_numeric(df.loc[df['Gene']==gene]['MAX_AF'],errors='coerce').replace('NaN',0).min()
        if max_af * min_af > 0.02:
          df = df.loc[df['Gene']!=gene]

      df.to_csv("$code.AR_comp.txt",index=False,sep='\t')

    workbook = Workbook("$code.filtered.xlsx")
    workbook.use_zip64()

    for file in os.listdir('.'):
      if file.startswith("$code") and file.endswith(".txt"):
        wb = workbook
        tab = csv.reader(open(''+file, 'rt'), delimiter='\t')
        ws = wb.add_worksheet(file.split('.')[1])
        for row, data in enumerate(tab):
          ws.write_row(row,0,data)

    workbook.close()
    EOF
  >>>

  runtime {
    docker: "mawumag/anno-vep"
    cpu: "1"
    memory: "3 GiB"
  }
  output {
    File annotated_tsv = code + ".anno.txt.zip"
    File filtered_xlsx = code + ".filtered.xlsx"
  }
}