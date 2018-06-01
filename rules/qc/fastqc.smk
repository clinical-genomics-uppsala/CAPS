# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq

rule fastqc:
    input:
        lambda wildcards: get_fastq(wildcards,units,'fq1'),
        lambda wildcards: get_fastq(wildcards,units,'fq2')
    output:
        html="fastqc/{sample}-{unit}.html",
        zip="fastqc/{sample}-{unit}.zip"
    log:
        "logs/fastqc/{sample}-{unit}.log"
    threads: 1
    wrapper:
        "0.19.3/bio/fastqc"
