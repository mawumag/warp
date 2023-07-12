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

import "PairedFastQsToUnmappedBam.wdl" as FastqToBam
import "../broad/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" as ExomeGermlineSingleSample

workflow FastQsToGvcf {
  input {
    String sample_id
    Array[File]+ fastq_1
    Array[File]+ fastq_2
    Array[Pair[File,File]] fastq = zip(fastq_1, fastq_2)
    PapiSettings papi_settings
    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    File target_interval_list
    File bait_interval_list
    String bait_set_name
  }

  scatter(fastq_pair in fastq) {
    call FastqToBam.ConvertPairedFastQsToUnmappedBamWf {
      input:
        sample_name = sample_id,
        fastq_1 = fastq_pair.left,
        fastq_2 = fastq_pair.right,
        readgroup_name = basename(fastq_pair.left),
        library_name = sample_id,
        platform_unit = "NA",
        run_date = "1970-01-01",
        platform_name = "ILLUMINA",
        sequencing_center = "NA"
    }
  }

  call ExomeGermlineSingleSample.ExomeGermlineSingleSample {
    input:
      papi_settings = papi_settings,
      references = references,
      scatter_settings = scatter_settings,
      target_interval_list = target_interval_list,
      bait_interval_list = bait_interval_list,
      bait_set_name = bait_set_name,
      sample_and_unmapped_bams = {
        "sample_name": sample_id,
        "base_file_name": sample_id,
        "flowcell_unmapped_bams": ConvertPairedFastQsToUnmappedBamWf.output_unmapped_bam,
        "final_gvcf_base_name":  sample_id,
        "unmapped_bam_suffix": ".unmapped.bam"
      }
    
  }

  output {
    File aligned_cram = ExomeGermlineSingleSample.output_cram
    File aligned_cram_index = ExomeGermlineSingleSample.output_cram_index
    File output_vcf = ExomeGermlineSingleSample.output_vcf
  }
}