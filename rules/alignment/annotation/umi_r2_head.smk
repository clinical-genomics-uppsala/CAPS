# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#rule annotate_bam_with_umis:
#    input:
#        bam="mapped/{sample}.{unit}.primerclip.bam",
#        umi="umi/{sample}.{unit}.R2.UMIs.fastq"
#    output:
#        bam="mapped/{sample}.{unit}.primerclip.umi.bam"
#    log:
#        "logs/umis/{sample}.{unit}.fgbioAnnoBam.txt"
#    wrapper:
#        "fgbio/bio/fgbio/annotatebamwithumis"

rule annotate_bam_with_umis_part:
    input:
        bam="mapped/{sample}.{unit}.{part}.primerclip.bam",
        umi="umi/{sample}.{unit}.{part}.R2.UMIs.fastq"
    output:
        bam="mapped/{sample}.{unit}.{part}.primerclip.umi.bam"
    log:
        "logs/umis/{sample}.{unit}.{part}.fgbioAnnoBam.txt"
    wrapper:
        "fgbio/bio/fgbio/annotatebamwithumis"
