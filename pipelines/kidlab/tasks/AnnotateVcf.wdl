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

task EnsemblVepAnnotateVcf {
  input {
    File input_vcf
    String output_vcf_basename
    File cadd_snv
    File cadd_snv_index
    File cadd_indel
    File cadd_indel_index
    File gnomad_exomes
    File gnomad_exomes_index
    File gerp_scores
    File alphamissense
    File alphamissense_index
    File domino
    File msc
    File gdi
    File connectome
    File PID_panel
    File PID_extra
    Int cpu = 2
    String memory = "8 GiB"
  }

  command <<<
    set -o pipefail
    set -e
    
    vep \
      --offline \
      -i ~{input_vcf} \
      -o STDOUT \
      --per_gene \
      --everything \
      --plugin CADD,~{cadd_snv},~{cadd_indel} \
      --plugin Conservation,~{gerp_scores} \
      --plugin AlphaMissense,file=~{alphamissense} \
      --custom ~{gnomad_exomes},gnomAD,vcf,exact,0,nhomalt \
      --force_overwrite \
      --fork ~{cpu} \
      --vcf | \
    bcftools +anno-vep - -- Domino "~{domino}" | \
    bcftools +anno-vep - -- MSC "~{msc}" | \
    bcftools +anno-vep - -- GDI "~{gdi}" | \
    bcftools +anno-vep - -- Connectome "~{connectome}" | \
    bcftools +anno-vep - -- PID "~{PID_panel}" | \
    bcftools +anno-vep -Oz -o ~{output_vcf_basename}.anno.vcf.gz - -- PIDg "~{PID_extra}"
  >>>

  runtime {
    docker: "mawumag/anno-vep-db"
    cpu: cpu
    memory: memory
    disks: "local-disk 300 SSD"
    bootDiskSizeGb: 50
  }
  output {
    File output_vcf = "~{output_vcf_basename}.anno.vcf.gz"
  }
}