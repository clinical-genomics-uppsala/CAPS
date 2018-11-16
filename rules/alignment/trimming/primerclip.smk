# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#rule sort_preprimerclip_queryname:
#    input:
#        "mapped/{sample}.{unit}.bam"
#    output:
#        temp("mapped/{sample}.{unit}.qsorted.sam")
#    log:
#        "logs/umi/qsort/{sample}.{unit}.log"
#    params:
#        sort_order="queryname",
#        extra="VALIDATION_STRINGENCY=LENIENT"
#    wrapper:
#        "bio/picard/sortsam"

rule sort_preprimerclip_queryname:
    input:
        "mapped/{sample}.{unit}.{part}.bam"
    output:
        temp("mapped/{sample}.{unit}.{part}.qsorted.sam")
    log:
        "logs/umi/qsort/{sample}.{unit}.{part}.log"
    params:
        sort_order="queryname",
        extra="VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "master/bio/picard/sortsam"

#rule primerclip:
#    input:
#        sam="mapped/{sample}.{unit}.qsorted.sam",
#        master_file="master.bed"
#    output:
#        sam=temp("mapped/{sample}.{unit}.primerclip.sam")
#    log:
#        "logs/mapped/primerclip/{sample}.{unit}.log"
#    wrapper:
#        "bio/primerclip"

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

#rule primerclip_bam_generation:
#    input:
#        "mapped/{sample}.{unit}.primerclip.sam"
#    output:
#      "mapped/{sample}.{unit}.primerclip.bam"
#    params:
#        "-Sb" # optional params string
#    wrapper:
#        "bio/samtools/view"

rule primerclip_bam_generation:
    input:
        "mapped/{sample}.{unit}.{part}.primerclip.sam"
    output:
      temp("mapped/{sample}.{unit}.{part}.primerclip.bam")
    params:
        "-Sb" # optional params string
    wrapper:
        "master/bio/samtools/view"

#rule sort_preUMI_queryname_part:
#    input:
#        "mapped/{sample}.{unit}.{part}.bam"
#    output:
#        temp("mapped/{sample}.{unit}.{part}.qsorted.sam")
#    log:
#        "logs/umi/qsort/{sample}.{unit}.{part}.log"
#    params:
#        sort_order="queryname",
#        extra="VALIDATION_STRINGENCY=LENIENT"
#    wrapper:
#        "bio/picard/sortsam"
