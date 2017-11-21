# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule queryname_sort_reads:
  input:
      "mapped/{sample}.sorted.bam"
  output:
      bam = "mapped/{sample}.queryname_sorted.bam"
  wrapper:
      "https://raw.githubusercontent.com/clinical-genomics-uppsala/snakemake-wrappers/master/bio/sort/queryname/wrapper.py"

rule amplicon_mapping:
  input:
      "mapped/{sample}.queryname_sorted.bam"
  output:
      bam = "mapped/{sample}.amplicon_annotated.bam",
      bed = "amplicon_information/{sample}.amplicon_annotated.bed"
  log:
      "logs/amplicon_annotated/{sample}.amplicon_annotated.log"
  params:
      path_gatk = config.get('path_gatk',"gatk.jar"),
      genome_ref = config['reference_genome'],
      design_file = "/DiagnosticPanel_Lung_20160222.selection.bed"
  wrapper:
      "https://raw.githubusercontent.com/clinical-genomics-uppsala/snakemake-wrappers/fix-logging/bio/amplicon_mapping/wrapper.py"

rule coordinate_sort_amplicon_mapped_reads:
    input:
        "mapped/{sample}.amplicon_annotated.bam"
    output:
        bam = "mapped/{sample}.amplicon_annotated.sorted.bam"
    threads: 3
    wrapper:
        "0.19.3/bio/samtools/sort"

rule create_bam_index_for_amplicon_mapped_bam:
    input:
        "mapped/{sample}.amplicon_annotated.sorted.bam"
    output:
        "mapped/{sample}.amplicon_annotated.sorted.bam.bai"
    params:
        ""
    wrapper:
        "0.17.4/bio/samtools/index"
