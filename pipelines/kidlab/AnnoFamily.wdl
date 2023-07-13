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
import "tasks/AnnotateVcf.wdl" as AnnotateVcf
import "tasks/GeneratePedFile.wdl" as GeneratePedFile
import "tasks/FilterVcf.wdl" as FilterVcf

workflow AnnoFamily {
  input {
    String family_id
    Array[String]+ samples
    Array[String]+ kinship
    Array[String]+ sex
    File input_vcf
    File vep_archive
    String vep_data_dir = "vep_data"
    File cadd_snv
    File cadd_snv_index
    File cadd_indel
    File cadd_indel_index
    File gnomad_exomes
    File gnomad_exomes_index
    File gerp_scores
  }

  call SplitVcf.SplitVcf {
    input:
      family_id = family_id,
      samples = samples,
      kinship = kinship,
      input_vcf = input_vcf
  }
  
  call AnnotateVcf.EnsemblVepAnnotateVcf {
    input:
      input_vcf = SplitVcf.output_vcf,
      output_vcf_basename = family_id,
      vep_archive = vep_archive,
      vep_data_dir = vep_data_dir,
      cadd_snv = cadd_snv,
   	  cadd_snv_index = cadd_snv_index,
      cadd_indel = cadd_indel,
      cadd_indel_index = cadd_indel_index,
      gnomad_exomes = gnomad_exomes,
      gnomad_exomes_index = gnomad_exomes_index,
      gerp_scores = gerp_scores
  }

  call GeneratePedFile.GeneratePedFile {
    input:
      family_id = family_id,
      samples = samples,
      kinship = kinship,
      sex = sex
  }

  call FilterVcf.FilterVcf {
    input:
      input_vcf = EnsemblVepAnnotateVcf.output_vcf,
      ped_file = GeneratePedFile.ped_file,
      code = family_id
  }

  output {
    File annotated_vcf = EnsemblVepAnnotateVcf.output_vcf
    File annotated_tsv = FilterVcf.annotated_tsv
    File filtered_xlsx = FilterVcf.filtered_xlsx
  }
}