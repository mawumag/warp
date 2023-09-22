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
  }

  call AnnotateVcf.EnsemblVepAnnotateVcf {
    input:
      input_vcf = input_vcf,
      output_vcf_basename = cohort_id,
  }

  output {
    File annotated_vcf = EnsemblVepAnnotateVcf.output_vcf
  }
}