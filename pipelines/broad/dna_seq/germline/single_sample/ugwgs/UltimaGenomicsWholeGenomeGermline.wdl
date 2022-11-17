version 1.0

import "../../../../../../pipelines/broad/dna_seq/somatic/single_sample/ugwgs/UltimaGenomicsWholeGenomeCramOnly.wdl" as UltimaGenomicsWholeGenomeCramOnly
import "../../../../../../tasks/broad/UltimaGenomicsWholeGenomeGermlineTasks.wdl" as Tasks
import "../../../../../../tasks/broad/Utilities.wdl" as Utilities
import "../../../../../../tasks/broad/GermlineVariantDiscovery.wdl" as VariantDiscoverTasks
import "../../../../../../tasks/broad/UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates.wdl" as UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates
import "../../../../../../tasks/broad/InternalTasks.wdl" as InternalTasks
import "../../../../../../tasks/broad/Qc.wdl" as QC
import "../../../../../../tasks/broad/UltimaGenomicsWholeGenomeGermlineQC.wdl" as UltimaGenomicsWholeGenomeGermlineQC
import "../../../../../../structs/dna_seq/UltimaGenomicsWholeGenomeGermlineStructs.wdl" as Structs
import "../../../../../../pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl" as ReblockGVCF


workflow UltimaGenomicsWholeGenomeGermline {
  input {

    ContaminationSites contamination_sites
    AlignmentReferences alignment_references
    VariantCallingSettings variant_calling_settings
    VcfPostProcessing vcf_post_processing

    # Sample Information
    Array[File]? input_cram_list
    Array[File]? input_bam_list
    String base_file_name

    Float rsq_threshold = 1.0
    Boolean merge_bam_file = true
    Boolean make_haplotype_bam = false
    Int reads_per_split = 20000000
    String filtering_model_no_gt_name = "rf_model_ignore_gt_incl_hpol_runs"

    File original_gvcf
    File original_gvcf_index
    String flow_order
  }

  meta {
    allowNestedInputs: true
  }

  parameter_meta {
    contamination_sites: "Struct containing files for contamination estimation"
    alignment_references: "Struct containing reference files for alignment with BWA mem"
    variant_calling_settings: "Struct containing reference files for variant calling with HaplotypeCaller"
    vcf_post_processing: "Struct containing reference files for VCF post-processing: annotation and filtering"
    input_cram_list: "Array of CRAM files to be used as workflow input. Must be specified if `input_bam_list` is not provided"
    input_bam_list: "Array of unmapped BAM files to be used as workflow input. Must be specified if `input_cram_list` is not provided"
    base_file_name: "Base name for each of the output files."
    rsq_threshold: "Threshold for a read quality metric that is produced by the sequencing platform"
    merge_bam_file: "Boolean indicating if by-interval bamout files from HaplotypeCaller should be merged into a single BAM"
    reads_per_split: "Number of reads by which to split the CRAM prior to alignment"
    filtering_model_no_gt_name: "String describing the optional filtering model; default set to rf_model_ignore_gt_incl_hpol_runs"
  }

  String pipeline_version = "1.0.6"


  References references = alignment_references.references

  call RemoveTreeScoreAnnotation {
    input:
      original_gvcf = original_gvcf,
      original_gvcf_index = original_gvcf_index
  }

  call Tasks.ConvertGVCFtoVCF {
    input:
      input_gvcf         = RemoveTreeScoreAnnotation.output_vcf,
      input_gvcf_index   = RemoveTreeScoreAnnotation.output_vcf_index,
      output_vcf_name    = base_file_name + '.vcf.gz',
      references         = references
  }

  # VCF post-processings
  call Tasks.AnnotateVCF {
    input :
      input_vcf               = ConvertGVCFtoVCF.output_vcf,
      input_vcf_index         = ConvertGVCFtoVCF.output_vcf_index,
      references              = references,
      reference_dbsnp         = vcf_post_processing.ref_dbsnp,
      reference_dbsnp_index   = vcf_post_processing.ref_dbsnp_index,
      flow_order              = flow_order,
      final_vcf_base_name     = base_file_name
  }

  call Tasks.AddIntervalAnnotationsToVCF {
    input:
      input_vcf             = AnnotateVCF.output_vcf_annotated,
      input_vcf_index       = AnnotateVCF.output_vcf_annotated_index,
      final_vcf_base_name   = base_file_name,
      annotation_intervals  = vcf_post_processing.annotation_intervals
  }

  call Tasks.TrainModel {
    input:
      input_file                = AddIntervalAnnotationsToVCF.output_vcf,
      input_file_index          = AddIntervalAnnotationsToVCF.output_vcf_index,
      input_vcf_name            = base_file_name,
      blocklist_file            = vcf_post_processing.training_blocklist_file,
      ref_fasta                 = references.ref_fasta,
      ref_index                 = references.ref_fasta_index,
      runs_file                 = vcf_post_processing.runs_file,
      apply_model               = filtering_model_no_gt_name,
      annotation_intervals      = vcf_post_processing.annotation_intervals,
      exome_weight              = vcf_post_processing.exome_weight,
      exome_weight_annotation   = vcf_post_processing.exome_weight_annotation
  }


  call Tasks.AnnotateVCF_AF {
    input :
      input_vcf             = AddIntervalAnnotationsToVCF.output_vcf,
      input_vcf_index       = AddIntervalAnnotationsToVCF.output_vcf_index,
      af_only_gnomad        = vcf_post_processing.af_only_gnomad,
      af_only_gnomad_index  = vcf_post_processing.af_only_gnomad_index,
      final_vcf_base_name   = base_file_name
  }

  call Tasks.FilterVCF {
    input:
      input_vcf               = AnnotateVCF_AF.output_vcf_annotated,
      input_model             = select_first([vcf_post_processing.filtering_model_no_gt,TrainModel.model_pkl]),
      runs_file               = vcf_post_processing.runs_file,
      references              = references,
      model_name              = filtering_model_no_gt_name,
      filter_cg_insertions    = vcf_post_processing.filter_cg_insertions,
      final_vcf_base_name     = base_file_name,
      flow_order              = flow_order,
      annotation_intervals    = vcf_post_processing.annotation_intervals
  }

  call Tasks.MoveAnnotationsToGvcf {
    input:
      filtered_vcf        = FilterVCF.output_vcf_filtered,
      filtered_vcf_index  = FilterVCF.output_vcf_filtered_index,
      gvcf                = original_gvcf,
      gvcf_index          = original_gvcf_index
  }

  call ReblockGVCF.ReblockGVCF {
    input:
      gvcf = MoveAnnotationsToGvcf.output_gvcf,
      gvcf_index = MoveAnnotationsToGvcf.output_gvcf_index,
      calling_interval_list = variant_calling_settings.wgs_calling_interval_list,
      ref_dict = alignment_references.references.ref_dict,
      ref_fasta = alignment_references.references.ref_fasta,
      ref_fasta_index = alignment_references.references.ref_fasta_index,
      tree_score_cutoff = vcf_post_processing.remove_low_tree_score_sites_cutoff,
      annotations_to_keep_command = vcf_post_processing.annotations_to_keep_command_for_reblocking
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_gvcf = ReblockGVCF.output_vcf
    File output_gvcf_index = ReblockGVCF.output_vcf_index
    File output_vcf = ConvertGVCFtoVCF.output_vcf
    File output_vcf_index = ConvertGVCFtoVCF.output_vcf_index

    # VCF post-processing
    File filtered_vcf = FilterVCF.output_vcf_filtered
    File filtered_vcf_index = FilterVCF.output_vcf_filtered_index
  }

}


task RemoveTreeScoreAnnotation {
  input {
    File original_gvcf
    File original_gvcf_index
    String base_file_name
  }

  command {
    set -eo pipefail

    bcftools annotate -x INFO/TREE_SCORE -o ~{base_file_name}.g.vcf.gz -O z ~{original_gvcf}
    bcftools index -t ~{base_file_name}.g.vcf.gz
  }
  output {
    File output_vcf = "~{base_file_name}.g.vcf.gz"
    File output_vcf_index = "~{base_file_name}.g.vcf.gz.tbi"
  }

  runtime {
        memory: "8GB"
        disks: "local-disk " + (ceil(size(original_gvcf, "GB")) * 2 + 10) + " HDD"
        docker: "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.6-1.10.2-0.1.16-1659548257"
        cpu: 1
    }
}