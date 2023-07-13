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

task SplitVcf {
  input {
    String family_id
    Array[String]+ samples
    Array[String]+ kinship
    File input_vcf
  }


  command <<<
    set -o pipefail
    set -e

    echo -e ~{sep="," samples} "\t" ~{sep="," kinship} "\t" ~{family_id} > sample.list
    bcftools +split ~{input_vcf} -i'GT="alt"' -S sample.list -o . -Oz
  >>>

  runtime {
    docker: "mawumag/anno-vep"
    cpu: "1"
    memory: "3 GiB"
  }

  output {
    File output_vcf = family_id + ".vcf.gz"
  }
}