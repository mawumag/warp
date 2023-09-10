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

import "tasks/SplitVcf.wdl" as SplitVcf
import "tasks/GeneratePedFile.wdl" as GeneratePedFile
import "tasks/CalculateGenotypePosteriors.wdl" as CalculateGenotypePosteriors
import "tasks/FilterVcf.wdl" as FilterVcf

workflow FilterFamily {
  input {
    String family_id
    Array[String]+ samples
    Array[String]+ kinship
    Array[String]+ sex
    File annotated_cohort_vcf
    Array[File] supporting_callsets
    Array[File] supporting_callsets_indices
  }

  call SplitVcf.SplitVcf {
    input:
      family_id = family_id,
      samples = samples,
      kinship = kinship,
      input_vcf = annotated_cohort_vcf
  }
  
  call GeneratePedFile.GeneratePedFile {
    input:
      family_id = family_id,
      samples = samples,
      kinship = kinship,
      sex = sex
  }

  call CalculateGenotypePosteriors.CalculateGenotypePosteriors {
    input :
      input_vcf = SplitVcf.output_vcf,
      ped_file = GeneratePedFile.ped_file,
      supporting_callsets = supporting_callsets,
      supporting_callsets_indices = supporting_callsets_indices,
      output_vcf_basename = family_id
  }

  call FilterVcf.FilterVcf {
    input:
      input_vcf = CalculateGenotypePosteriors.output_vcf,
      ped_file = GeneratePedFile.ped_file,
      code = family_id
  }

  output {
    File annotated_vcf = CalculateGenotypePosteriors.output_vcf
    File annotated_tsv = FilterVcf.annotated_tsv
    File filtered_xlsx = FilterVcf.filtered_xlsx
  }
}