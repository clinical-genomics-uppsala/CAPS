# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""

Expects the global variable config of at least the following structure
...............................................................................
---
illuminaclip_file: /path/to/illumina.fa
cutadapt:
  adapters_1:
    first_pair_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC,
    second_pair_adapter: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
...............................................................................

Expects a samples config
With the following dictionary structure (panel, sample_number and lane are
required)
samples = {
    'cutadapt': {sample_id: value},
    'fq1': {samplei_id: value},
    'fq2' {sample_id: value}
}

Example of a sample.tsv that can be imported using pandas (columns need to be
tab separated)
...............................................................................
sample     adapters     fq1                fq2
sample1    adapters_1   sample1.R1.fastq   sample1.R2.fastq
sample2    adapters_1   sample2.R1.fastq   sample2.R2.fastq
...............................................................................
"""
from scripts.lib.common.utils import get_fastq, sample_id, reverse_complement

rule cutadapt:
   input:
       [lambda wildcards: get_fastq(wildcards,units,'fq1'),
        lambda wildcards: get_fastq(wildcards,units,'fq2')]
   output:
       fastq1="trimmed/{sample}.{unit}.R1.cutadapt.fastq.gz",
       fastq2="trimmed/{sample}.{unit}.R2.cutadapt.fastq.gz",
       qc = "logs/trimmed/{sample}.{unit}.cutadapt.qc.txt"
   params:
       " --minimum-length 1",
       lambda wildcards: " -a " + config["cutadapt"][samples['adapters'][sample_id(wildcards)]]["first_pair_adapter"],
       lambda wildcards: " -A " + reverse_complement(config["cutadapt"][samples['adapters'][sample_id(wildcards)]]["second_pair_adapter"])
   wrapper:
       "0.17.4/bio/cutadapt/pe"
