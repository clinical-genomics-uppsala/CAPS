# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

_amplicon_mapping_input = "mapped/{sample}.sorted.bam"
try:
    _amplicon_mapping_input = amplicon_mapping_input
except:
    pass

_amplicon_mapping_output = "mapped/{sample}.amplicon_annotated.sorted.bam"
try:
    _amplicon_mapping_output = amplicon_mapping_output
except:
      pass

rule queryname_sort_reads:
  input:
      _amplicon_mapping_input
  output:
      bam = temp("mapped/{sample}.queryname_sorted.bam")
  params:
      remove_secondary_alignment = True
  wrapper:
      "master/bio/sort/queryname"

rule amplicon_mapping:
  input:
      "mapped/{sample}.queryname_sorted.bam"
  output:
      bam = temp("mapped/{sample}.amplicon_annotated.bam"),
      bed = "amplicon_information/{sample}.amplicon_annotated.bed"
  params:
      path_gatk = config.get('path_gatk',"gatk.jar"),
      genome_ref = config['reference_genome'],
      design_file =lambda wildcards: samples['amplicon_file'][wildcards.sample]
  wrapper:
      "master/bio/amplicon_mapping"

rule coordinate_sort_amplicon_mapped_reads:
    input:
        "mapped/{sample}.amplicon_annotated.bam"
    output:
        _amplicon_mapping_output
    wrapper:
        "0.19.3/bio/samtools/sort"

rule create_bam_index_for_amplicon_mapped_bam:
    input:
        _amplicon_mapping_output
    output:
        _amplicon_mapping_output + ".bai"
    params:
        ""
    wrapper:
        "0.27.1/bio/samtools/index"
