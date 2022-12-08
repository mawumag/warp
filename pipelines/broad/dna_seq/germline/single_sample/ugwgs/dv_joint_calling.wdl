version 1.0

workflow dv_joint_calling {
    input {
        Array[File] gvcf_gz
        File bed_file
        String output_name
    }

    call GLnexus { input:
        gvcf = gvcf_gz,
        bed_file = bed_file,
        config = "DeepVariant",
        output_name = output_name
    }

   output {
        File pvcf_gz = GLnexus.pvcf_gz
    }
}

task GLnexus {
    input {
        Array[File]+ gvcf
        File bed_file
        String config
        String output_name
    }

    command <<<
        set -ex -o pipefail
        export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.1
        numactl --interleave=all glnexus_cli --config "~{config}" --list --bed "~{bed_file}" "~{write_lines(gvcf)}" | bcftools view - | bgzip -@ 4 -c > "~{output_name}.vcf.gz"
    >>>

    runtime {
        docker: "quay.io/mlin/glnexus:v0.6.0-0-ge7bbd0c"
        disks: "local-disk 400 HDD"
        cpu: "16"
    }

    output {
        File pvcf_gz = "${output_name}.vcf.gz"
    }
}