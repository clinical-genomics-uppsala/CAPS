# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

_import_read_input = lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
try:
    _import_read_input = import_read_input
except:
    pass

_import_read_output = "trimmed/{sample}.{unit}.{read}.fastq.gz"
try:
    _import_read_output = import_read_output
except:
    pass

rule import_read:
    input:
        _import_read_input
    output:
        _import_read_output
    shell:
        "cp {input} {output}"
