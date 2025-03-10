# 1.2.24
2023-06-23 (Date of Last Commit)

* Added a move command to work around a maximum file path length in featureCounts
* Updated STARsolo version to v2.7.10b for the StarsoloFastq task; does not impact this workflow
* Updated STARsolo argument for counting mode to GeneFull_Ex50pAS; does not impact this workflow

# 1.2.23 
2023-05-04 (Date of Last Commit)

* Updated the CheckInputs WDL for the Optimus workflow. This changes does impact snSS2

# 1.2.22

2023-04-23 (Date of Last Commit)

* Updated the STARalign task; does not affect this workflow

# 1.2.21
2023-04-19 (Date of Last Commit)

* Updated warp-tools docker which included a fix for a small bug in create_snrna_optimus.py that was causing the script not to run


# 1.2.20
2023-03-27 (Date of Last Commit)

* SlideSeq-specific and Optimus-specific changes to Metrics.wdl. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline

# 1.2.19
2023-03-15 (Date of Last Commit)

* SlideSeq-specific and Optimus-specific changes to Metrics.wdl. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline
* Updated warp-tools docker to support the Optimus and SlideSeq changes 

# 1.2.18
2023-02-28 (Date of Last Commit)

* Added a new task to the workflow that reads the tar_star_reference file to obtain the genomic reference source, build version, and annotation version and outputs the information as a txt file

# 1.2.17
2023-02-13 (Date of Last Commit)

* SlideSeq-specific changes to LoomUtils.wdl. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline

# 1.2.16
2023-02-07 (Date of Last Commit)

* Updated the input checks for the Optimus pipeline task; this has no effect on this pipeline
* Added disk to the checkInputArrays task 

# 1.2.15
2023-01-23 (Date of Last Commit)

* Added 'Disk' to task runtime sections to support running on Azure
* Addressed mb/gb memory specification inconsistencies in LoomUtils and CheckInput

# 1.2.14
2022-09-20 (Date of Last Commit)

* Updated tasks in StarAlign.wdl to use an updated STAR docker image.

# 1.2.13

2022-09-01 (Date of Last Commit)

* Updated CheckInputs.wdl to use a lightweight alpine-bash image.

# 1.2.12
2022-08-31 (Date of Last Commit)

* Updated CountAlignments to use an updated docker image.

# 1.2.11
2022-08-23 (Date of Last Commit)

* Removed an unused script in pytools docker image.

# 1.2.10
2022-08-16 (Date of Last Commit)

* Updated LoomUtils.wdl to use a consolidated python utilities docker image. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline.

# 1.2.9
2022-08-08 (Date of Last Commit)

* Updated TrimAdapters runtime docker URL.

# 1.2.8
2022-07-21 (Date of Last Commit)

* Updated STARsoloFastq runtime docker URL.

# 1.2.7
2022-05-18 (Date of Last Commit)

* Updated merge npz docker in StarAlign.wdl

# 1.2.6
2022-04-22 (Date of Last Commit)

* Updated LoomUtils.wdl for a task in the Optimus pipeline. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline.

# 1.2.5
2022-04-19 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.2.4

2022-04-06 (Date of Last Commit)

* Updated STARsoloFastq task in StarAlign.wdl; this change does not affect the Mutl-snSS2 workflow.

# 1.2.3
2022-02-25 (Date of Last Commit)

* Updated LoomUtils.wdl for a task in the Optimus pipeline. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline.

# 1.2.2
2022-02-10 (Date of Last Commit)

* Rebuilt a docker to merge outputs of STAR in StarAlign.wdl task and moved it to a public location.

# 1.2.1
2022-02-07 (Date of Last Commit)

* Updated a task in STARalign.wdl related to the Optimus pipeline. This pipeline has not been changed.

# 1.2.0
2022-01-07 (Date of Last Commit)

* Fixed missing metadata issue in the loom file
# 1.1.2
2021-11-19 (Date of Last Commit)

* Updated STARsoloFastq to use 'output_bam_basename' to name the aligned bam. This is consistent with versions 4.2.7 and older. This change has no impact MultiSampleSmartSeq2SingleNucleus
# 1.1.1
2021-11-15 (Date of Last Commit)
* Updated remove-reads-on-junctions.py in the FeatureCounts.wdl to use python3 instead of python2.

# 1.1.0
2021-09-16 (Date of Last Commit)

* Removed the Smart-seq2 Single Nucleus workflow (SmartSeq2SingleNucleus.wdl) and changed the workflow tasks to run multiple samples in the same VM. This change is expected to make the pipeline faster and cheaper.
* Renamed the StarAlignFastq.StarAlignFastqPairedEnd task to StarAlign.StarAlignFastqMultisample

# 1.0.4
2021-09-10 (Date of Last Commit)

Added the option "--soloBarcodeReadLength 0" STARsoloFastq task to support alignment in Optimus. This change has no impact on MultiSampleSmartSeq2SingleNucleus

# 1.0.3
2021-09-02 (Date of Last Commit)

* Added a new StarSolo task for Optimus in the StarAlign.wdl. However, the same wdl
  contains other Star tasks that are used in the smartseq2 single nuclei for paired and 
  single stranged fastq files. As a result, the smartseq2 processing is not expected to 
  change. 

# 1.0.2
2021-08-02 (Date of Last Commit)

* Increased the version number to make new release tag for Dockstore 

# 1.0.1
2021-07-19 (Date of Last Commit)

* Updated SmartSeq2 to accommodate spaces in input_name

# 1.0.0

2021-05-17 (Date of First Commit)

* This is the first release of the Smart-seq2 Multi-Sample Single Nuclei workflow
