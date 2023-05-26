version 1.0

workflow BamMetrics {
  input {
    File input_bam
    File annotation_gtf
    String output_name
    String mt_sequence
  }
  call DropseqMetrics {
	  input:
      input_bam = input_bam,
      annotation_gtf = annotation_gtf,
      output_name = output_name,
      mt_sequence = mt_sequence
  }
  
  output {
    File final_out = DropseqMetrics.metric_output
  }

}
task DropseqMetrics {
  input {
    File input_bam
    File annotation_gtf
    String output_name
    String mt_sequence
  }
  command <<<
    /usr/gitc/Drop-seq_tools-2.5.3/SingleCellRnaSeqMetricsCollector \
    --ANNOTATIONS_FILE ~{annotation_gtf} \
    --INPUT ~{input_bam} \
    --NUM_CORE_BARCODES 40000 \
    --OUTPUT ~{output_name} \
    --MT_SEQUENCE ~{mt_sequence} \
    --CELL_BARCODE_TAG CB \
    --TMP_DIR /data/docker_tmp/tmp/
  >>>
  runtime {
	docker: "ekiernan/samtools-picard:v1"
  }
  output {
    File metric_output = "~{output_name}"
  }
}