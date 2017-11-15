# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule queryname_sort_reads:
  input:
    "mapped/{sample}.bam"
  output:
    bam = "amplicon_mapped/{sample}.queryname_sorted.bam"
  wrapper:
    "https://raw.githubusercontent.com/clinical-genomics-uppsala/snakemake-wrappers/master/bio/sort/queryname/wrapper.py"

rule amplicon_mapping:
  input:
    "amplicon_mapped/{sample}.queryname_sorted.bam"
  output:
    bam = "amplicon_mapped/{sample}.amplicon_annotated.bam",
    bed = "amplicon_mapped/{sample}.amplicon_annotated.bed"
  params:
    path_gatk = "/GenomeAnalysisTKLite_molecules.jar",
    genome_ref = "/hg19.with.mt.fasta",
    design_file = "/DiagnosticPanel_Lung_20160222.selection.bed"
  wrapper:
    "https://raw.githubusercontent.com/clinical-genomics-uppsala/snakemake-wrappers/master/bio/amplicon_mapping/wrapper.py"

rule coordinate_sort_amplicon_mapped_reads:
    input:
        "amplicon_mapped/{sample}.amplicon_annotated.bam"
    output:
        bam = "amplicon_mapped/{sample}.sorted.bam"
    wrapper:
        "https://raw.githubusercontent.com/clinical-genomics-uppsala/snakemake-wrappers/master/bio/sort/coordinate/wrapper.py"

rule create_bam_index_for_amplicon_mapped_bam:
    input:
        "amplicon_mapped/{sample}.sored.bam"
    output:
        "amplicon_mapped/{sample}.sorted.bam.bai"
    params:
        ""
    wrapper:
        "0.17.4/bio/samtools/index"
