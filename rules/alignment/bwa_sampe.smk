# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.common.utils import get_fastq_files, get_now

rule bwa_alignment:
    input:
        [lambda wildcards: "trimmed/" + get_fastq_files(wildcards,samples,"fq1"),
         lambda wildcards: "trimmed/" + get_fastq_files(wildcards,samples,"fq2")]
    output:
        temp("mapped/{sample}-{unit}.unsorted.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
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
        "mapped/{sample}-{unit}.unsorted.bam"
    output:
        bam = "mapped/{sample}-{unit}.sorted.bam"
    wrapper:
        "https://raw.githubusercontent.com/clinical-genomics-uppsala/snakemake-wrappers/master/bio/sort/coordinate/wrapper.py"

rule create_bam_index:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        "mapped/{sample}-{unit}.sorted.bam.bai"
    params:
        ""
    wrapper:
        "0.17.4/bio/samtools/index"
