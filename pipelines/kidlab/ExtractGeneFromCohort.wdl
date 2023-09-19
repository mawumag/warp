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

workflow ExtractGeneFromCohort {
  input {}

  call ExtractGene {}

  output {
    File gene_tsv = ExtractGene.gene_tsv
  }
}

task ExtractGene {
  input {
    File annotated_cohort_vcf
    String gene_name
    String gene_id
  }

  command <<<
    input=~{annotated_cohort_vcf}
    output="~{gene_name}.txt"
    gene_id=~{gene_id}
    format='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t[%SAMPLE %GT, ]\t'
    csq=$(bcftools +split-vep -l "$input" | cut -f 2 | tr '\n' '\t' | sed 's/\t$/\\n/' | sed 's/\t/\\t%/g' | sed 's/^/%/')
    format+=$csq

    bcftools +split-vep -g <(echo "${gene_id}") -H -f"$format" -i'GT="alt"' $input | \
    sed -r '1 s/\[[0-9]+\]//g' | \
    sed -r '1 s/^#[ ]?//' > $output
  >>>

  runtime {
    docker: "mawumag/anno-vep"
    cpu: "1"
    memory: "3 GiB"
  }

  output {
    File gene_tsv = gene_name + ".txt"
  }
}