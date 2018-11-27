# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq_files, get_now

_bwa_alignment_input = lambda wildcards: ["trimmed/" + get_fastq_files(wildcards,samples,'fq1'), "trimmed/" + get_fastq_files(wildcards,samples,'fq2')]
try:
    if bwa_mem_input:
        _bwa_alignment_input = lambda wildcards: bwa_mem_input
except:
      pass

_bwa_alignment_output = "mapped/{sample}.{unit}.{part}.bam"
try:
  _bwa_alignment_output = bwa_mem_output
except:
  pass

rule bwa_alignment:
    input:
        reads=_bwa_alignment_input
    output:
        _bwa_alignment_output
    log:
        "logs/bwa_mem/{sample}.{unit}.{part}.log"
    threads: 3
    params:
        index=config['reference_genome'],
        extra=lambda wildcards: r"-M -R '@RG\tID:" + get_now() + "_" + wildcards.sample + r"\tSM:" + wildcards.sample + r"\tPL:illumina'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-@ 16"
    wrapper:
        "0.27.1/bio/bwa/mem"
