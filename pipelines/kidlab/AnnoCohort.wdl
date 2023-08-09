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

import "tasks/AnnotateVcf.wdl" as AnnotateVcf

workflow AnnoCohort {
  input {
    String cohort_id
    File input_vcf
    File cadd_snv
    File cadd_snv_index
    File cadd_indel
    File cadd_indel_index
    File gnomad_exomes
    File gnomad_exomes_index
    File gerp_scores
    File domino
    File msc
    File gdi
    File connectome
    File PID_panel
    File PID_extra
  }

  call AnnotateVcf.EnsemblVepAnnotateVcf {
    input:
      input_vcf = input_vcf,
      output_vcf_basename = cohort_id,
      cadd_snv = cadd_snv,
   	  cadd_snv_index = cadd_snv_index,
      cadd_indel = cadd_indel,
      cadd_indel_index = cadd_indel_index,
      gnomad_exomes = gnomad_exomes,
      gnomad_exomes_index = gnomad_exomes_index,
      gerp_scores = gerp_scores,
      domino = domino,
      msc = msc,
      gdi = gdi,
      connectome = connectome,
      PID_panel = PID_panel,
      PID_extra = PID_extra
  }

  output {
    File annotated_vcf = EnsemblVepAnnotateVcf.output_vcf
  }
}