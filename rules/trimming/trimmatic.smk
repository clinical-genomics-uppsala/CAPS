# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

_head_crop_input = lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
try:
  _head_crop_input = head_crop_input
except:
  pass

_head_crop_output = temp("trimmed/.temp/{sample}.{unit}.{read}.fastq")
try:
  _head_crop_output = head_crop_output
except:
  pass

rule head_crop:
    input:
        _head_crop_input
    output:
        _head_crop_output
    run:
        if input[0].endswith("gz"):
        shell("zcat {input} > {output}")
     else:
        shell("cat {input} > {output}")
