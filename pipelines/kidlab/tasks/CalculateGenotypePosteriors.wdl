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

task CalculateGenotypePosteriors {
  input {
    File input_vcf
    File ped_file
    Array[File] supporting_callsets
    Array[File] supporting_callsets_indices
    String output_vcf_basename
  }

  command <<<
    gatk IndexFeatureFile -I ~{input_vcf}
    gatk --java-options "-Xms3G -Xmx4G" CalculateGenotypePosteriors \
      -V ~{input_vcf} \
      --supporting-callsets ~{sep="," supporting_callsets} \
      --ped ~{ped_file} \
      -O ~{output_vcf_basename}.vcf.gz
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    cpu: 1
    memory: "5 GiB"
  }

  output {
    File output_vcf = output_vcf_basename + ".vcf.gz"
  }
}