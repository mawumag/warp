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
#   Tasks used in multiple workflows
# CHANGELOG
#

import "structs.wdl" as Structs

task ExtractSampleNameFlowOrder{
    input{
        File input_bam
        File monitoring_script
        Int preemptible_tries
        String docker
        References references
        Boolean no_address
        String gitc_path = "/usr/gitc"
    }
    parameter_meta {
        input_bam: {
            localization_optional: true
        }
    }

    command <<<
        set -e
        set -o pipefail

        bash ~{monitoring_script} > monitoring.log &

        java -jar ~{gitc_path}/GATK_ultima.jar  GetSampleName  \
            -I ~{input_bam} \
            -R ~{references.ref_fasta} \
            -O sample_name.txt

        gsutil cat ~{input_bam} | java -jar ~{gitc_path}/GATK_ultima.jar ViewSam -I /dev/stdin \
                --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/FO:/) {print substr($i,4,4)}}}' | head -1 \
                > flow_order.txt

        gsutil cat ~{input_bam} | java -jar ~{gitc_path}/GATK_ultima.jar ViewSam -I /dev/stdin \
                --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/BC:/) {print substr($i,4)}}}' | head -1 \
                > barcode.txt

        gsutil cat ~{input_bam} | java -jar ~{gitc_path}/GATK_ultima.jar ViewSam -I /dev/stdin \
                --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/ID:/) {print substr($i,4)}}}' | head -1 \
                > id.txt

    >>>

    runtime {
        cpu: 1
        memory: "2 GB"
        preemptible: preemptible_tries
        noAddress: no_address
        docker: docker
    }

    output {
        String sample_name = read_string("sample_name.txt")
        String flow_order = read_string("flow_order.txt")
        String barcode_seq = read_string("barcode.txt")
        String readgroup_id = read_string("id.txt")
        File sample_name_file = "sample_name.txt"
        File flow_order_file = "flow_order.txt"
        File barcode_seq_file = "barcode.txt"
        File readgroup_id_file = "id.txt"
        File monitoring_log = "monitoring.log"
    }
}

# Aggregate picard metrics and coverage and variant calling reports into one h5
task AggregateMetrics {
    input {
        File monitoring_script
        Array[File]? picard_files
        String base_file_name
        File? coverage_all_h5
        File? short_report_h5
        File? extended_report_h5
        File? no_gt_report_h5
        Float? contamination_stdout
        String docker
        Int disk_size = 20
        Int preemptible_tries
        Boolean no_address
    }

    command <<<
    set -eo pipefail
    bash ~{monitoring_script} > monitoring.log &
    source ~/.bashrc
    conda activate genomics.py3

    collect_existing_metrics.py \
        ~{true="--metric_files " false="" defined(picard_files)} ~{ sep=' ' picard_files} \
        ~{"--coverage_h5 " + coverage_all_h5} \
        ~{"--short_report_h5 " + short_report_h5} \
        ~{"--extended_report_h5 " + extended_report_h5} \
        ~{"--no_gt_report_h5 " + no_gt_report_h5} \
        ~{"--contamination_stdout " + contamination_stdout} \
        --output_h5 ~{base_file_name}.aggregated_metrics.h5
    >>>

    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File monitoring_log = "monitoring.log"
        File aggregated_metrics_h5 = "~{base_file_name}.aggregated_metrics.h5"
    }
}

# Convert the aggregated metrics h5 file to json
task ConvertAggregatedMetricsToJson {
    input {
        File monitoring_script
        File? aggregated_metrics_h5
        String base_file_name
        String docker
        Int preemptible_tries
        Boolean no_address
        Int disk_size = 20
    }

    command <<<
    set -eo pipefail

    bash ~{monitoring_script} > monitoring.log &
    source ~/.bashrc
    conda activate genomics.py3

    convert_h5_to_json.py \
            --root_element "metrics" \
            --ignored_h5_key_substring histogram \
            --input_h5 ~{aggregated_metrics_h5} \
            --output_json ~{base_file_name}.aggregated_metrics.json

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "1 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File aggregated_metrics_json = "~{base_file_name}.aggregated_metrics.json"
    }
}

task ScatterIntervalList {
    input {
        File monitoring_script
        File interval_list
        Int scatter_count
        Int break_bands_at_multiples_of
        String docker
        String gitc_path
        Boolean no_address
        String dummy_input_for_call_caching  # !UnusedDeclaration
    }
    command <<<
    bash ~{monitoring_script} > monitoring.log &
    echo ~{dummy_input_for_call_caching}
    set -eo pipefail
    mkdir out
    java -Xms4g -jar ~{gitc_path}picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{break_bands_at_multiples_of} \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i+1).rjust(len(str(len(intervals))),"0") + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
    >>>
    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int(stdout())
        File monitoring_log = "monitoring.log"
    }
    runtime {
        memory: "6 GB"
        docker: docker
        noAddress: no_address
    }
}

task IntervalListOfGenome {
  input  {
    File ref_fai
    File ref_dict
    String docker
    Int disk_size
    Int preemptible_tries
    Boolean no_address
    File monitoring_script
  }
  command <<<
    set -e
    bash ~{monitoring_script} > monitoring.log &

    cat ~{ref_fai} | grep -v '[ME_-]' | awk '{print $1"\t"1"\t"$2"\t+\t."}' > modifai
    cat ~{ref_dict} modifai > genome.interval_list
  >>>
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    continueOnReturnCode: true
    maxRetries: 1
    noAddress: no_address
  }
  output{
    File monitoring_log = "monitoring.log"
    File interval_list = "genome.interval_list"
  }
}

task ConcatMetricsJsons {
    input {
        File monitoring_script
        Array[File] jsons
        String base_file_name
        String docker
        Int preemptible_tries
        Boolean no_address
        Int disk_size = 20
    }

    command <<<
    set -eo pipefail

    bash ~{monitoring_script} > monitoring.log &
    source ~/.bashrc
    conda activate genomics.py3

    echo ~{sep=',' jsons} > json_files.txt

    python <<CODE
    import json

    json_out = {'metrics': {}}
    with open('json_files.txt') as file_handler:
        json_files = file_handler.read().rstrip().split(',')
    for jf in json_files:
        with open(jf) as file_handler:
            json_content = json.load(file_handler)
        if "metrics" not in json_content:
            raise ValueError(f"Invalid json file {jf} - does not contain a 'metrics' key")
        joint_keys = set(json_content["metrics"]).intersection(set(json_out["metrics"]))
        if len(joint_keys) > 0:
            raise ValueError(f"Encountered duplicate metrics - {joint_keys}")
        json_out['metrics'].update(json_content['metrics'])
    with open('~{base_file_name}.aggregated_metrics.json','w') as file_handler:
        json.dump(json_out, file_handler, indent=2)
    CODE

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "1 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File aggregated_metrics_json = "~{base_file_name}.aggregated_metrics.json"
    }
}

task DownsampleCramBam {
    input {
        String base_file_name
        File input_cram_bam
        Float downsample_frac
        Int seed = 789
        References references
        File monitoring_script
        String docker
        Boolean no_address
        Int preemptibles
        Int disk_size = ceil((1.1+downsample_frac)*size(input_cram_bam,"GB") + 20)
    }
    String output_cram_name = base_file_name + ".cram"
    command <<<
        set -eo pipefail

        bash ~{monitoring_script} > monitoring.log &

        source ~/.bashrc
        conda activate genomics.py3

        samtools view \
          -s ~{seed}~{downsample_frac} \
          -C \
          -T ~{references.ref_fasta} \
          -o ~{output_cram_name} \
          --threads 4 \
          ~{input_cram_bam}

        samtools index ~{output_cram_name}
    >>>
    runtime {
        disks: "local-disk " + disk_size + " HDD"
        cpu: 4
        memory: "1 GB"
        preemptible: preemptibles
        docker: docker
        noAddress: no_address
    }
    output {
        File output_cram = "~{output_cram_name}"
        File output_cram_index = "~{output_cram_name}.crai"
        File monitoring_log = "monitoring.log"
    }
}

task ToFastq {
    input {
        File input_cram
        String base_file_name
        References references
        File monitoring_script
        String docker
        Boolean no_address
        Int preemptibles
        Int disk_size = ceil(3*size(input_cram,"GB") + 20)
    }
    String output_fq_name = base_file_name + ".fq.gz"
    command <<<
        set -eo pipefail

        bash ~{monitoring_script} > monitoring.log &

        source ~/.bashrc
        conda activate genomics.py3

        samtools fastq \
          --reference ~{references.ref_fasta}\
          -F 0x900 \
          -0 ~{output_fq_name} \
          --threads 2 \
          ~{input_cram}
    >>>
    runtime {
        disks: "local-disk " + disk_size + " HDD"
        cpu: 2
        memory: "1 GB"
        preemptible: preemptibles
        docker: docker
        noAddress: no_address
    }
    output {
        File output_fastq = "~{output_fq_name}"
        File monitoring_log = "monitoring.log"
    }
}

task ConcatHtmls {
   input {
        # input for concat
        Array[File] htmls

        # call arguments
        File monitoring_script
        String base_file_name

        # runtime arguments
        Int preemptible_tries
        String docker
        Int disk_size = ceil(2*size(htmls,"GB") + 1)



    }
    command <<<
        set -eo pipefail

        start=$(date +%s)
        source ~/.bashrc
        conda activate genomics.py3

        bash ~{monitoring_script} > monitoring.log &

        echo ~{sep=',' htmls} > htmls_files.txt
        ls -lstr
        cat htmls_files.txt

        i=1
        i=1; for file in $(cat htmls_files.txt| sed "s/,/ /g"); do cp $file $i.html ;
        cat $i.html >> ~{base_file_name}_aggregated.html; rm $i.html ;i=$((i + 1)) ; done

        ls -lstr
        echo "**************** D O N E ****************"

        end=$(date +%s)
        mins_elapsed=$(( ($end - $start) / 60))
        secs_elapsed=$(( ($end - $start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
            secs_elapsed=0$secs_elapsed
        fi
        echo "Total run time: $mins_elapsed:$secs_elapsed"
        ls -ltrsa

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        docker: docker
        cpu: "1"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        noAddress: true
    }

    output {
        File report_html = "~{base_file_name}_aggregated.html"
        File monitoring_log = "monitoring.log"
    }
}

