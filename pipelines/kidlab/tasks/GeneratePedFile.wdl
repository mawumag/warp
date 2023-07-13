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

task GeneratePedFile {
  input {
    String family_id
    Array[String]+ samples
    Array[String]+ kinship
    Array[String]+ sex
  }

  command <<<
    samples=(~{sep=" " samples})
    kinship=(~{sep=" " kinship})
    sex=(~{sep=" " sex})

    father="0"
    mother="0"

    arrayLength=${#samples[@]}

    for (( i=0; i<$arrayLength; i++ )); do
      if [[ ${kinship[$i]} == "father" ]]; then
        father=${samples[$i]}
      elif [[ ${kinship[$i]} == "mother" ]]; then
        mother=${samples[$i]}
      fi
    done

    for (( i=0; i<$arrayLength; i++ )); do
      if [[ ${sex[$i]} == "M" ]]; then
        sex_individual="1"
      elif [[ ${sex[$i]} == "F" ]]; then
        sex_individual="2"
      else
        sex_individual="3"
      fi

      if [[ ${kinship[$i]} == "father" ]]; then
        echo "~{family_id} ${samples[$i]} 0 0 1 1" >> ~{family_id}.ped
      elif [[ ${kinship[$i]} == "mother" ]]; then
        echo "~{family_id} ${samples[$i]} 0 0 2 1" >> ~{family_id}.ped
      elif [[ ${kinship[$i]} == "proband" ]]; then
        echo "~{family_id} ${samples[$i]} $father $mother $sex_individual 2" >> ~{family_id}.ped
      else
        echo "~{family_id} ${samples[$i]} $father $mother $sex_individual 1" >> ~{family_id}.ped
      fi
    done
  >>>

  runtime {
    docker: "mawumag/anno-vep"
    cpu: "1"
    memory: "3 GiB"
  }

  output {
    File ped_file = family_id + ".ped"
  }
}