# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule import_read:
    input:
        lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
    output:
        "trimmed/{sample}.{unit}.{read}.fastq.gz"
    shell:
        "cp {input} {output}"
