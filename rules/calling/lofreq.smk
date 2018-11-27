# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


_lofreq_input = "mapped/consensus/{sample}.consensus.bam"
try:
    _lofreq_input = lofreq_input
except:
    pass

_lofreq_output = "mapped/consensus/{sample}.consensus.bam"
try:
    _lofreq_output = lofreq_output
except:
    pass


rule lofreq:
    input:
        _lofreq_input
    output:
        temp("calls/{sample}.lofreq.vcf")
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
        _lofreq_output
    wrapper:
        "lofreq2indelovlp-wrapper/bio/lofreq/tools/lofreq2indelovlp"
