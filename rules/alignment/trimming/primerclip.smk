# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

_primerclip_input = lambda wildcards: ["trimmed/" + get_fastq_files(wildcards,samples,'fq1'), "trimmed/" + get_fastq_files(wildcards,samples,'fq2')]
try:
    _primerclip_input = primerclip_input
except:
    pass

_primerclip_output = "mapped/{sample}.{unit}.{part}.primerclip.bam"
try:
    _primerclip_output = primerclip_output
  except:
    pass

rule sort_preprimerclip_queryname:
    input:
        _primerclip_input
    output:
        temp("mapped/{sample}.{unit}.{part}.qsorted.sam")
    log:
        "logs/umi/qsort/{sample}.{unit}.{part}.log"
    params:
        sort_order="queryname",
        extra="VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "master/bio/picard/sortsam"

rule primerclip:
    input:
        sam="mapped/{sample}.{unit}.{part}.qsorted.sam",
        master_file="master.bed"
    output:
        sam=temp("mapped/{sample}.{unit}.{part}.primerclip.sam")
    log:
        "logs/mapped/primerclip/{sample}.{unit}.{part}.log"
    wrapper:
        "primerclip/bio/primerclip"

rule primerclip_bam_generation:
    input:
        "mapped/{sample}.{unit}.{part}.primerclip.sam"
    output:
        _primerclip_output
    params:
        "-Sb" # optional params string
    wrapper:
        "master/bio/samtools/view"
