version 1.0

task global {
  command {}
  runtime { docker: "ubuntu:focal" }
  output {
    String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698"
    String gitc_jar_path = "/usr/gitc/"
    String ug_vc_docker = "gcr.io/ultima-data-307918/ugvc:ugvc_50512d"
    String ug_gatk_picard_docker = "gcr.io/ganymede-331016/ug_gatk_picard:0.8"
    String dv_docker = "gcr.io/terra-project-249020/deepvariant:ug-1.4.4_8ebb16b5"
    String vbid_docker = "us.gcr.io/broad-gotc-prod/verify-bam-id:f6cb51761861e57c43879aa262df5cf8e670cf7c-1606775309"
    String sentieon_docker = "gcr.io/ganymede-331016/sentieon_docker:sentieon-genomics-202112.06"
    String ua_docker = "gcr.io/ultima-data-307918/ua:ua_394e0a"
    String crammer_docker = "gcr.io/ganymede-331016/crammer:1.3_0c527b2e"
    String trimmer_docker = "gcr.io/ganymede-331016/trimmer:vc_1.2.0"
    String single_cell_docker = "gcr.io/ganymede-331016/single_cell_10x:gatk-v0.7rc3_scp-v1.12"
    String starsolo_docker = "gcr.io/ganymede-331016/star:2.7.8a"
    String fastqc_docker = "quay.io/biocontainers/fastqc:0.11.9--0"
    String bwa_meth_docker = "gcr.io/ganymede-331016/nfcore_methylseq:1.6.1"
    String monitoring_script = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
  }
}
