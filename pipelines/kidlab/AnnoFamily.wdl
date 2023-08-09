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
import "tasks/CalculateGenotypePosteriors.wdl" as CalculateGenotypePosteriors
import "tasks/FilterVcf.wdl" as FilterVcf

workflow AnnoFamily {
  input {
    String family_id
    Array[String]+ samples
    Array[String]+ kinship
    Array[String]+ sex
    File unannotated_cohort_vcf 
    File? annotated_cohort_vcf
    Array[File] supporting_callsets
    Array[File] supporting_callsets_indices
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

  call SplitVcf.SplitVcf {
    input:
      family_id = family_id,
      samples = samples,
      kinship = kinship,
      input_vcf = select_first([annotated_cohort_vcf, unannotated_cohort_vcf])
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

  if (!defined(annotated_cohort_vcf)) {
    call AnnotateVcf.EnsemblVepAnnotateVcf {
      input:
        input_vcf = CalculateGenotypePosteriors.output_vcf,
        output_vcf_basename = family_id,
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
  }

  call FilterVcf.FilterVcf {
    input:
      input_vcf = select_first([EnsemblVepAnnotateVcf.output_vcf, CalculateGenotypePosteriors.output_vcf]),
      ped_file = GeneratePedFile.ped_file,
      code = family_id
  }

  output {
    File annotated_vcf = select_first([EnsemblVepAnnotateVcf.output_vcf, CalculateGenotypePosteriors.output_vcf])
    File annotated_tsv = FilterVcf.annotated_tsv
    File filtered_xlsx = FilterVcf.filtered_xlsx
  }
}