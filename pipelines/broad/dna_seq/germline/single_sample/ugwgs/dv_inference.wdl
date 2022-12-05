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
# Runs DeepVariant calling pipeline for Ultima Genomics data.

import "../../../../../../tasks/broad/structs.wdl" as Structs
import "../../../../../../tasks/broad/general_tasks.wdl" as UGGeneralTasks
import "../../../../../../tasks/broad/globals.wdl" as Globals

workflow DVInference {
  input {
    # Workflow args
    String pipeline_version = "1.2.0"
    String base_file_name

    # Mandatory inputs
    Array[File] bam_files
    Array[File] bam_index_files
    References references

    # Scatter interval list args
    Int num_shards
    Int scatter_intervals_break # Maximal resolution for scattering intervals

    # Make examples args
    Float min_fraction_hmer_indels
    Float min_fraction_non_hmer_indels
    Float min_fraction_snps
    Int min_base_quality
    Int pileup_min_mapping_quality
    Int candidate_min_mapping_quality
    Int dbg_min_base_quality
    Int min_windows_distance
    String alt_aligned_pileup
    String? ug_channels_args
    File? subset_interval_list
    Int vsc_max_background_count = -1
    Float max_background_fraction = -1
    Boolean make_gvcf
    Float p_error
    # Call variants args
    File model
    File model_index
    File model_meta

    Int preemptible_tries = 1

    # Systematic error correction args
    File? background_bam_file
    File? background_bam_index_file

    String dummy_input_for_call_caching = ""

   # winval validations
   #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
   #@wv defined(background_bam_file) <-> defined(background_bam_index_file)
   #@wv prefix(model_index) == prefix(model_meta)
   #@wv prefix(model_index) == prefix (model)
   #@wv min_fraction_hmer_indels <= 1 and min_fraction_hmer_indels >= 0
   #@wv min_fraction_non_hmer_indels <= 1 and min_fraction_non_hmer_indels >= 0
   #@wv len(suffix(bam_files) & {".bam", ".cram"}) > 0
   #@wv len(suffix(bam_index_files) & {".bai", ".crai"}) > 0
   #@wv prefix(bam_index_files) == bam_files
   #@wv len(bam_files) >= 0
   #@wv len(bam_index_files) >= 0
   #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
   #@wv suffix(references['ref_dict']) == '.dict'
   #@wv suffix(references['ref_fasta_index']) == '.fai'
   #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
  }

  call Globals.global

  if (!defined(subset_interval_list)){
    call UGGeneralTasks.IntervalListOfGenome as IntervalListOfGenome{
      input:
        ref_fai = references.ref_fasta_index,
        ref_dict = references.ref_dict,
        disk_size = 1,
        preemptible_tries = preemptible_tries,
        docker = global.gitc_docker,
        monitoring_script = global.monitoring_script,
        no_address = true
    }
  }

  File interval_list = select_first([IntervalListOfGenome.interval_list, subset_interval_list])

  String make_examples_executable = (if defined(background_bam_file) then "multisample_make_examples" else "make_examples")
  String output_prefix = base_file_name #basename(bam_files[0], ".cram")

  call UGGeneralTasks.ScatterIntervalList as ScatterIntervalList{
    input:
      interval_list = interval_list,
      scatter_count = num_shards,
      break_bands_at_multiples_of = scatter_intervals_break,
      dummy_input_for_call_caching = dummy_input_for_call_caching,
      docker = global.gitc_docker,
      gitc_path = global.gitc_jar_path,
      no_address = true,
      monitoring_script = global.monitoring_script
  }

  scatter (interval in ScatterIntervalList.out){
    call MakeInferenceExamples {
      input:
        interval = interval,
        total_number_of_shards = ScatterIntervalList.interval_count,
        references = references,
        bam_files = bam_files,
        bam_index_files = bam_index_files,
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
        make_gvcf = make_gvcf,
        p_error   = p_error,
        dv_docker = global.dv_docker,
        make_examples_executable = make_examples_executable,
        monitoring_script = global.monitoring_script,
        preemptible_tries = preemptible_tries
    }
  }

  call CallVariants {
    input:
      examples = MakeInferenceExamples.output_examples,
      model = model,
      model_index = model_index,
      model_meta = model_meta,
      dv_docker = global.dv_docker,
      call_variants_exec = "call_variants",
      total_number_of_shards = ScatterIntervalList.interval_count,
      monitoring_script = global.monitoring_script
  }

  if (make_gvcf){
    Array[File] gvcf_records = MakeInferenceExamples.gvcf_record
  }

  call PostProcessing {
    input:
      called_records = CallVariants.output_records,
      gvcf_records = gvcf_records,
      make_gvcf    = make_gvcf,
      total_number_of_shards = ScatterIntervalList.interval_count,
      ref = references.ref_fasta,
      ref_index = references.ref_fasta_index,
      dv_docker = global.dv_docker,
      output_prefix = output_prefix,
      monitoring_script = global.monitoring_script
  }

  # this is done to fix the caching issue
  if (make_gvcf){ 
    File gvcf_maybe = PostProcessing.gvcf_file
    File gvcf_index_maybe = PostProcessing.gvcf_file_index
  }

  output {
    File monitoring_log     = PostProcessing.monitoring_log
    File vcf_file           = PostProcessing.vcf_file
    File vcf_index          = PostProcessing.vcf_index
    File? output_gvcf       = gvcf_maybe
    File? output_gvcf_index = gvcf_index_maybe

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
    Boolean make_gvcf 
    Float? p_error
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
  String gvcf_filename = basename(interval, ".interval_list") + ".gvcf.tfrecord.gz"
  String gvcf_string = if make_gvcf then ("--p_error " + p_error + " --gvcf " + gvcf_filename)  else ""
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
      ~{gvcf_string} 
    echo "Finished make_examples task!"
    touch ~{gvcf_filename}
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
    File gvcf_record = gvcf_filename
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
    Int disk_size = ceil(2 * size(examples, 'GB') + 2 * size(model, 'GB'))
  }
  command <<<
    set -eo pipefail

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
    gpuType: "nvidia-tesla-p100"
    gpuCount: 1
    nvidiaDriverVersion: "515.65.01"
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
    Int total_number_of_shards
    Array[File]? gvcf_records

    Int disk_size = ceil(48 * size(called_records, "GB") +
                         size(ref, "GB") + 
                         (if make_gvcf then size(select_first([gvcf_records]), "GB") else 0) + 
                         4 + (if make_gvcf then 5 else 0))
    Int memory = ceil(64 * ( size(called_records, "GB")) + 
      (if make_gvcf then 10*size(select_first([gvcf_records]), "GB") else 0)) + 4

  }

  
  command <<<
      bash ~{monitoring_script} > monitoring.log &
      set -eo pipefail

      echo ~{sep="," gvcf_records} > files.txt
      # gvcf names should be of the form gvcf.tfrecord-{00001}-of-{00096}.gz, and not 
      # as we get them from MakeInferenceExamples. This process fixes the names
      # also if gvcf_records is empty - the previous echo will create an empty line that we skip
      python -c 'import os; \
      file_list=open("files.txt").read().rstrip().split(","); \
      nfiles_string=str(len(file_list)).zfill(5); \
      [os.symlink(filename, f"gvcf.tfrecord-{str(i).zfill(5)}-of-{nfiles_string}.gz") \
      for i,filename in enumerate(file_list) if filename.strip()]'

      ls *.gz || true # ignore *.gz missing for make_gvcf=False


      /opt/deepvariant/bin/postprocess_variants \
      --ref ~{ref} \
      --infile ~{called_records} \
      --outfile ~{output_prefix}.vcf.gz \
      ~{if make_gvcf then "--nonvariant_site_tfrecord_path gvcf.tfrecord@" + total_number_of_shards + ".gz" else ""} \
      ~{if make_gvcf then "--gvcf_outfile " + output_prefix + ".g.vcf.gz" else ""} \
      --vcf_stats_report=False

      touch ~{output_prefix}.g.vcf.gz
      touch ~{output_prefix}.g.vcf.gz.tbi
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
    File gvcf_file = '~{output_prefix}.g.vcf.gz'
    File gvcf_file_index = '~{output_prefix}.g.vcf.gz.tbi'
  }
}
