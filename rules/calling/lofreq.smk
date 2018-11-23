# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule lofreq:
    input:
        "mapped/consensus/{sample}.consensus.bam"
    output:
        "calls/{sample}.lofreq.vcf"
    params:
        ref=config['reference_genome'],
        extra="-l merged_targets.bed"
    threads: 2
    wrapper:
        "lofreq-wrapper/bio/lofreq/call"

rule merge_lofreq_indels:
    input:
        "calls/{sample}.lofreq.vcf"
    output:
        "calls/{sample}.lofreq.merged.vcf"
    wrapper:
        "lofreq2indelovlp-wrapper/bio/lofreq/tools/lofreq2indelovlp"
