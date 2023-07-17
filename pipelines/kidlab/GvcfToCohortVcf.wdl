version 1.0

## Copyright Marco Baggio, KU Leuven, 2023
## Copyright Broad Institute, 2018 (for imported workflows)
## For questions: kuleuven@mawumag.com
## 
## This WDL implements the Broad Institute ExomeGermlineSingleSample pipeline starting
## from fastq files instead of unmapped bams
## 
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "GenerateSampleMap.wdl" as GenerateSampleMap
import "../broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl" as JointGenotyping

workflow GvcfToCohortVcf {
  input {}

  call GenerateSampleMap.GenerateSampleMap {}

  call JointGenotyping.JointGenotyping {
    input:
      sample_name_map = GenerateSampleMap.sample_map
  }

  output {
    File cohort_vcf = JointGenotyping.output_vcfs[0]
    File cohort_vcf_index = JointGenotyping.output_vcf_indices[0]
  }
}