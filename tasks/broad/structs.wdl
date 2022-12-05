version 1.0
# LICENSE
#   Copyright 2022 Ultima Genomics
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
# DESCRIPTION
#   Structures definitions used mainly in ug_variant_calling_pipeline.
# CHANGELOG in reverse chronological order
#
struct ContaminationSites {
  String contamination_sites_path
  File contamination_sites_vcf
  File contamination_sites_vcf_index
}

struct References {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
}

struct CramBamFiles {
  File input_cram_bam_list
  File input_cram_bam_index
}

# BWA alignment
struct AlignmentReferences {
  References references
  File ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
}

# BWA-METH alignment
struct BwaMethReferences {
  File bwa_meth_ref_idx
}

# UA alignment
struct UaReferences {
  File ua_index
  File ref_alt
  String ua_extra_args
  Boolean v_aware_alignment_flag
}

# UA-METH alignment
struct UaMethReferences {
  File index_c2t
  File index_g2a
}

struct VariantCallingSettings {
  File wgs_calling_interval_list
  Int break_bands_at_multiples_of
  Int haplotype_scatter_count
}

struct VcfPostProcessing {
  Array[File] annotation_intervals
  File? filtering_model_no_gt
  File af_only_gnomad
  File af_only_gnomad_index
  Boolean filter_cg_insertions
  File? filtering_blocklist_file
  File? training_blocklist_file
  Int? exome_weight
  String? exome_weight_annotation
  File? interval_list_override
  File runs_file
  String? filtering_model_with_gt_name_override
  Float max_duplication_in_reasonable_sample
  Float max_chimerism_in_reasonable_sample
  File ref_dbsnp
  File ref_dbsnp_index
  File wgs_coverage_interval_list
}

struct SimpleReadTrimmingParameters {
  String trimming_tool  # allowed values cutadapt or trimmer
  String adapter_5p
  String adapter_3p
  # the 5p and 3p adapters are mandatory because the python script in trimmer requires it, can be changed if needed
  Int? umi_length_3p
  Int? umi_length_5p
  Float max_error_rate_5p
  Float max_error_rate_3p
  Int min_overlap_5p
  Int min_overlap_3p
  Int? min_insert_length
  Int? max_insert_length
  # in practice only trimmer supports it - current implementation of cutadapt can only discard through ClipReads and
  # it's not the behavior we want

  String untrimmed_reads_action
  # in practice only trimmer supports it - either "" (do nothing), "filter" (mark in sam flag) or "discard", cutadapt
  # only applies a tag to trimmed reads - XF:Z:A
}

struct PositionalErrorRateParams {
  Boolean allow_softlclip
  Int? position_plot_step_size
  Int? max_positions_plot
  Int? edist_threshold
  Int? xscore_threshold
}

struct ReferenceDbSnp {
    File ref_dbsnp
    File ref_dbsnp_index
}

struct FeatureMapParams {
  Int scatter_count
  Int min_mapq
  Int snv_identical_bases
  Int snv_identical_bases_after
  Int min_score
  Int limit_score
  String extra_args
}

struct StarsoloBamParams {
  Boolean save_bam
  Boolean? sort_bam_override
  String? limit_sort_ram_override
  Boolean? include_unmapped_override
  Boolean? transcript_override
}

struct RuntimeParams {
  Int? memory_gb_override
  String? disk_type_override
  Float? disk_size_gb_override
  Int? cpu_num_override
  Boolean? no_address_override
  Int? preemptible_tries_override
  Int? max_retries_override
  String? gitc_path_override
}

struct Adapters10x {
  String adapter_5p
  String adapter_3p
  String adapter_middle
  Float adapter_min_error_rate_5p
  Int adapter_min_overlap_5p
  Float adapter_min_error_rate_middle
  Int adapter_min_overlap_middle
  Float adapter_min_error_rate_3p
  Int adapter_min_overlap_3p
}


struct StarsoloScores {
  Float? starsolo_score_del_open_override
  Float? starsolo_score_ins_open_override
  Float? starsolo_score_del_base_override
  Float? starsolo_score_ins_base_override

}