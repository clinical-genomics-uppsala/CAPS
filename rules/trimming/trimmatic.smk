# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule head_crop:
    input:
        lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
    output:
        temp("trimmed/.temp/{sample}.{unit}.{read}.fastq")
    run:
        if input[0].endswith("gz"):
        shell("zcat {input} > {output}")
     else:
        shell("cat {input} > {output}")
