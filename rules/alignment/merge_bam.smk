# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_split_and_unit_part_files

_merge_bam_input = lambda wildcards: ["mapped/{}.{}.bam".format(wildcards.sample,unit_part) for unit_part in get_split_and_unit_part_files(wildcards,units,config)]
try:
    _merge_bam_input = merge_bam_input
except:
    pass

_merge_bam_output = "mapped/{sample}.sorted.bam"
try:
    _merge_bam_output = merge_bam_output
except:
    pass

rule merge_bam_files:
    input:
        _merge_bam_input
    output:
        temp("mapped/{sample}.merged.bam")
    threads: 8
    wrapper:
        "0.27.1/bio/samtools/merge"

rule coordinate_sort_mapped_reads_merged:
    input:
        "mapped/{sample}.merged.bam" #"mapped1/{sample}.{unit}.unsorted.bam"
    output:
        bam = _merge_bam_output
    threads: 3
    wrapper:
        "0.27.1/bio/samtools/sort"

rule create_bam_index:
    input:
        _merge_bam_output
    output:
        _merge_bam_output + ".bai"
    params:
        ""
    wrapper:
        "0.27.1/bio/samtools/index"
