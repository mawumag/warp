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

    #Temp DV inputs
    File? background_bam_file
    File? background_bam_index_file
    Float min_fraction_hmer_indels = 0.12
    Float min_fraction_non_hmer_indels = 0.06
    Float min_fraction_snps = 0.12
    Int min_base_quality = 5
    Int pileup_min_mapping_quality =5
    Int candidate_min_mapping_quality = 5
    Int dbg_min_base_quality = 0
    Int min_windows_distance = 20
    String alt_aligned_pileup = "none"
    String? ug_channels_args = "--aux_fields_to_keep tp,t0 --skip_bq_channel --channels hmer_deletion_quality,hmer_insertion_quality,non_hmer_insertion_quality"
    Int vsc_max_background_count = -1
    Float max_background_fraction = -1
    String dv_docker = "gcr.io/terra-project-249020/deepvariant:ug-1.4.4_8ebb16b5"
    String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698"
    String gitc_jar_path = "/usr/gitc/"
    File model = "gs://jukebox-broad-data-share/deepvariant/model/germline/v1.1/model.ckpt-740000.data-00000-of-00001"
    File model_index = "gs://jukebox-broad-data-share/deepvariant/model/germline/v1.1/model.ckpt-740000.index"
    File model_meta = "gs://jukebox-broad-data-share/deepvariant/model/germline/v1.1/model.ckpt-740000.meta"
    Int preemptible_tries = 1
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

  String pipeline_version = "1.0.5"

  String make_examples_executable = (if defined(background_bam_file) then "multisample_make_examples" else "make_examples")
  String output_prefix = base_file_name 
  String monitoring_script="gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"

  References references = alignment_references.references

  call UltimaGenomicsWholeGenomeCramOnly.UltimaGenomicsWholeGenomeCramOnly {
    input:
      contamination_sites = contamination_sites,
      alignment_references = alignment_references,
      input_cram_list = input_cram_list,
      input_bam_list = input_bam_list,
      base_file_name = base_file_name,
      rsq_threshold = rsq_threshold,
      reads_per_split = reads_per_split,
      vcf_post_processing = vcf_post_processing,
      save_bam_file = true
  }

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call Utilities.ScatterIntervalList {
    input:
      interval_list               = variant_calling_settings.wgs_calling_interval_list,
      scatter_count               = variant_calling_settings.haplotype_scatter_count,
      break_bands_at_multiples_of = variant_calling_settings.break_bands_at_multiples_of
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by factor 2 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int hc_divisor = ScatterIntervalList.interval_count / 2

  # Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    call MakeInferenceExamples {
      input:
        interval = ScatterIntervalList.out[index],
        total_number_of_shards = ScatterIntervalList.interval_count,
        references = references,
        bam_files = [select_first([UltimaGenomicsWholeGenomeCramOnly.output_bam])],
        bam_index_files = [select_first([UltimaGenomicsWholeGenomeCramOnly.output_bam_index])],
        background_bam_file = background_bam_file,
        background_bam_index_file = background_bam_index_file,
        alt_aligned_pileup = alt_aligned_pileup,
        min_fraction_hmer_indels = min_fraction_hmer_indels,
        min_fraction_non_hmer_indels = min_fraction_non_hmer_indels,
        min_fraction_snps = min_fraction_snps,
        min_base_quality = min_base_quality,
        candidate_min_mapping_quality = candidate_min_mapping_quality,
        pileup_min_mapping_quality = pileup_min_mapping_quality,
        dbg_min_base_quality = dbg_min_base_quality,
        min_windows_distance = min_windows_distance,
        vsc_max_background_count = vsc_max_background_count,
        vsc_max_background_fraction = vsc_max_background_count,
        ug_channels_args = ug_channels_args,
        dv_docker = dv_docker,
        make_examples_executable = make_examples_executable,
        monitoring_script = monitoring_script,
        preemptible_tries = preemptible_tries
    }
  }

  call CallVariants {
    input:
      examples = MakeInferenceExamples.output_examples,
      model = model,
      model_index = model_index,
      model_meta = model_meta,
      dv_docker = dv_docker,
      call_variants_exec = "call_variants",
      total_number_of_shards = ScatterIntervalList.interval_count,
      monitoring_script = monitoring_script
  }
  call PostProcessing {
    input:
      called_records = CallVariants.output_records,
      ref = references.ref_fasta,
      ref_index = references.ref_fasta_index,
      dv_docker = dv_docker,
      output_prefix = output_prefix,
      monitoring_script = monitoring_script
  }
  

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = PostProcessing.output_vcf
    File output_vcf_index = PostProcessing.output_vcf_index

    File output_cram = UltimaGenomicsWholeGenomeCramOnly.output_cram
    File output_cram_index = UltimaGenomicsWholeGenomeCramOnly.output_cram_index
    File output_cram_md5 = UltimaGenomicsWholeGenomeCramOnly.output_cram_md5

    File selfSM = UltimaGenomicsWholeGenomeCramOnly.selfSM
    Float contamination = UltimaGenomicsWholeGenomeCramOnly.contamination

    # STATISTIC COLLECTION
    File quality_yield_metrics = UltimaGenomicsWholeGenomeCramOnly.quality_yield_metrics
    File wgs_metrics = UltimaGenomicsWholeGenomeCramOnly.wgs_metrics
    File raw_wgs_metrics = UltimaGenomicsWholeGenomeCramOnly.raw_wgs_metrics
    File duplicate_metrics = UltimaGenomicsWholeGenomeCramOnly.duplicate_metrics
    File agg_alignment_summary_metrics = UltimaGenomicsWholeGenomeCramOnly.agg_alignment_summary_metrics
    File? agg_alignment_summary_pdf = UltimaGenomicsWholeGenomeCramOnly.agg_alignment_summary_pdf
    File agg_gc_bias_detail_metrics = UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_summary_metrics
    File agg_quality_distribution_pdf = UltimaGenomicsWholeGenomeCramOnly.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = UltimaGenomicsWholeGenomeCramOnly.agg_quality_distribution_metrics
    Float duplication_rate = UltimaGenomicsWholeGenomeCramOnly.duplication_rate
    Float chimerism_rate = UltimaGenomicsWholeGenomeCramOnly.chimerism_rate
    Boolean is_outlier_data = UltimaGenomicsWholeGenomeCramOnly.is_outlier_data

    String sample_name = UltimaGenomicsWholeGenomeCramOnly.sample_name
    String flow_order = UltimaGenomicsWholeGenomeCramOnly.flow_order
    String barcode = UltimaGenomicsWholeGenomeCramOnly.barcode
    String id = UltimaGenomicsWholeGenomeCramOnly.id
  }

}

task MakeInferenceExamples {
  input {
    Array[File] bam_files
    Array[File] bam_index_files
    File interval
    Float total_number_of_shards
    References references

    String alt_aligned_pileup
    Float min_fraction_hmer_indels
    Float min_fraction_non_hmer_indels
    Float min_fraction_snps
    Int min_base_quality
    Int candidate_min_mapping_quality
    Int pileup_min_mapping_quality
    Int dbg_min_base_quality
    Int min_windows_distance
    String? ug_channels_args
    Int vsc_max_background_count
    Float vsc_max_background_fraction

    # Background sample inputs
    File? background_bam_file
    File? background_bam_index_file
    String make_examples_executable

    String dv_docker
    File monitoring_script

    Int preemptible_tries
  }
  # Estimate output_size that fits candidate generated parameters (assuming constant image size)
  #   More sensitive thresholds (such as used for somatic variant detection) yield more examples (images)
  #   which consume more disk-space, regardless of the input-size.
  Float min_threshold_tmp = if min_fraction_hmer_indels < min_fraction_non_hmer_indels then min_fraction_hmer_indels else min_fraction_non_hmer_indels
  Float min_threshold = if min_threshold_tmp < min_fraction_snps then min_threshold_tmp else min_fraction_snps
  # solving linear equation for min_threshold=0.03 -> output_size=3000GB,
  #                             min_threshold=0.12 -> output_size=60GB
  Float a = -98000 / 3
  Float b = 3980
  Float expected_genome_wide_output = min_threshold * a + b
  Int expected_output_size = ceil(expected_genome_wide_output / total_number_of_shards)
  Int inputs_size =  ceil((size(bam_files, "GB") + size(background_bam_file, "GB")) / total_number_of_shards + size(references.ref_fasta, "GB"))
  
  Float c_i = 1.2  # inputs safety factor
  Float c_o = 2.5 # outputs safety factor
  Int disk_size = ceil(c_i * inputs_size + c_o * expected_output_size)

  String examples_filename = basename(interval, ".interval_list") + ".inference_examples.tfrecord.gz"

  parameter_meta {
      bam_files: {
          localization_optional: true
      }
      background_bam_file: {
          localization_optional: true
      }
  }

  command <<<
    set -eo pipefail

    bash ~{monitoring_script} > monitoring.log &

    /opt/gatk-4.2.6.1/gatk --java-options "-Xms1G" PrintReads \
        -I ~{sep=' -I ' bam_files} \
        -O input.bam \
        -L ~{interval} \
        -R ~{references.ref_fasta}

    defined_background=~{true="true" false="false" defined(background_bam_file)}
    if $defined_background ; then
      /opt/gatk-4.2.6.1/gatk --java-options "-Xms1G" PrintReads \
          -I ~{background_bam_file} \
          -O background.bam \
          -L ~{interval} \
          -R ~{references.ref_fasta}
    fi

    /opt/gatk-4.2.6.1/gatk IntervalListToBed -I ~{interval} -O interval.bed

    /opt/deepvariant/bin/~{make_examples_executable} \
      --mode calling \
      --ref ~{references.ref_fasta} \
      --reads "input.bam~{true=";background.bam" false="" defined(background_bam_file)}" \
      --examples ~{examples_filename} \
      --alt_aligned_pileup ~{alt_aligned_pileup} \
      --min_base_quality ~{min_base_quality} \
      --dbg_min_base_quality ~{dbg_min_base_quality} \
      --vsc_min_fraction_indels ~{min_fraction_non_hmer_indels} \
      --vsc_min_fraction_hmer_indels ~{min_fraction_hmer_indels} \
      --vsc_min_fraction_snps ~{min_fraction_snps} \
      --ws_min_windows_distance ~{min_windows_distance} \
      --vsc_max_background_count ~{vsc_max_background_count} \
      --vsc_max_background_fraction ~{vsc_max_background_fraction} \
      --min_mapping_quality ~{pileup_min_mapping_quality} \
      --candidate_min_mapping_quality ~{candidate_min_mapping_quality} \
      --output_only_sample_role_to_train \
      --regions interval.bed \
      ~{ug_channels_args} \

    echo "Finished make_examples task!"

  >>>
  runtime {
    memory: "2.5 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: dv_docker
    preemptible: preemptible_tries
  }
  output {
    File monitoring_log = "monitoring.log"
    File output_examples = examples_filename
  }
}

task CallVariants{
  input{
    Array[File] examples
    File model
    File model_index
    File model_meta
    String dv_docker
    String call_variants_exec
    Int total_number_of_shards
    File monitoring_script
    Int disk_size = ceil(1.2 * size(examples, 'GB') + 2 * size(model, 'GB'))
  }
  command <<<
    bash ~{monitoring_script} > monitoring.log &

    # Create symlinks with the naming convention prefered by deepvariant
    echo ~{sep="," examples} > files.txt
    python -c 'import os; \
      file_list=open("files.txt").read().rstrip().split(","); \
      nfiles_string=str(len(file_list)).zfill(5); \
      [os.symlink(filename, f"examples.tfrecord-{str(i).zfill(5)}-of-{nfiles_string}.gz") \
      for i,filename in enumerate(file_list)]'

    ls *.gz

    /opt/deepvariant/bin/~{call_variants_exec} \
      --outfile call_variants_output.tfrecord.gz \
      --examples examples.tfrecord@~{total_number_of_shards}.gz \
      --checkpoint $(echo ~{model_index} | sed s/.index$// )
  >>>
  runtime {
    memory: "14 GB"
    cpu: "4"
    disks: "local-disk " + disk_size + " HDD"
    docker: dv_docker
    gpuType: "nvidia-tesla-k80"
    gpuCount: 1
    nvidiaDriverVersion: "418.87.00"
  }
  output {
    File monitoring_log = "monitoring.log"
    File output_records = 'call_variants_output.tfrecord.gz'
  }
}

task PostProcessing{
  input{
    File called_records
    File ref
    File ref_index
    String dv_docker
    String output_prefix
    File monitoring_script
    Boolean make_gvcf = false
    Int disk_size = ceil(48 * size(called_records, "GB") +
                         size(ref, "GB") + 
                         (if make_gvcf then 12 else 0 ) + 4)
    Int memory = ceil(64 * size(called_records, "GB"))
  }
  command <<<
    bash ~{monitoring_script} > monitoring.log &
    /opt/deepvariant/bin/postprocess_variants \
      --ref ~{ref} \
      --infile ~{called_records} \
      --outfile ~{output_prefix}.vcf.gz \
      ~{true='--gvcf_outfile temp_name.g.vcf.gz' false='' make_gvcf} \
      --vcf_stats_report=False

    rename_gvcf=~{true="true" false="false" make_gvcf}
    if $rename_gvcf ; then
      mv temp_name.g.vcf.gz ~{output_prefix}.g.vcf.gz
    fi
  >>>
  runtime {
    memory: "~{memory} GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: dv_docker
  }
  output {
    File monitoring_log = "monitoring.log"
    File vcf_file = '~{output_prefix}.vcf.gz'
    File vcf_index = '~{output_prefix}.vcf.gz.tbi'
    File? gvcf_file = '~{output_prefix}.g.vcf.gz'
  }
}