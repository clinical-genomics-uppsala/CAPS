# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq_files, get_now


rule bwa_alignment:
    input:
        lambda wildcards: "trimmed/" + get_fastq_files(wildcards,samples,"fq1"),
        lambda wildcards: "trimmed/" + get_fastq_files(wildcards,samples,"fq2")
    output:
        temp("mapped/{sample}.{unit}.unsorted.bam")
    log:
        "logs/bwa_mem/{sample}.{unit}.log"
    threads: 8
    params:
        index=config['reference_genome'],
        extra=lambda wildcards: r"-M -R '@RG\tID:" + get_now() + "_" + wildcards.sample + r"\tSM:" + wildcards.sample + r"\tPL:illumina'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-@ 3"
    wrapper:
        "0.17.4/bio/bwa/mem"

rule coordinate_sort_mapped_reads:
    input:
        "mapped/{sample}.{unit}.unsorted.bam"
    output:
        bam = temp("mapped/{sample}.{unit}.sorted.bam")
    threads: 3
    wrapper:
        "0.19.3/bio/samtools/sort"

def get_units(wildcards,units):
    return [wildcards.sample +"-" + unit for unit in units.loc[wildcards.sample].index]

rule merge_bam_files:
    input:
        lambda wildcards: expand("mapped/{sample_units}.sorted.bam", sample_units=get_units(wildcards,units) )
    output:
        "mapped/{sample,[A-Za-z0-9-_]+}.sorted.bam"
    threads: 8
    wrapper:
        "0.19.3/bio/samtools/merge"

rule create_bam_index:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        "mapped/{sample,[A-Za-z0-9_-]+}.sorted.bam.bai"
    params:
        ""
    wrapper:
        "0.17.4/bio/samtools/index"
