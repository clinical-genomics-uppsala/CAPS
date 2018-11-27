# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq

_fastq_input = lambda wildcards: get_fastq(wildcards,units,'fq1'), lambda wildcards: get_fastq(wildcards,units,'fq2')
try:
    _fastq_input = fastq_input
except:
    pass

_fastq_html_output = "fastqc/{sample}.{unit}.html"
try:
    _fastq_html_output = fastq_html_output
except:
    pass

_fastq_zip_output = "fastqc/{sample}.{unit}.zip"
try:
    _fastq_zip_output = fastq_zip_output
except:
    pass

rule fastqc:
    input:
        _fastq_input
    output:
        html=_fastq_html_output,
        zip=_fastq_zip_output
    log:
        "logs/fastqc/{sample}.{unit}.log"
    wrapper:
        "0.27.1/bio/fastqc"
