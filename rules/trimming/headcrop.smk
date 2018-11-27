# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq

_headcrop_read_input = "trimmed/{sample}.{unit}.{part}.{read}.fastq.gz"
try:
    _headcrop_read_input = headcrop_read_input
except:
    pass

_headcrop_read_output = temp("trimmed/{sample}.{unit}.{part}.{read}.headcrop.fastq.gz")
try:
    _headcrop_read_output = headcrop_read_output
except:
    pass

rule headcrop_read:
   input:
        _headcrop_read_input
   output:
        _headcrop_read_output
   log:
      "logs/trimming/{sample}.{unit}.{part}.{read}.headcrop.trimmomatic.txt"
   params:
      trimmer=["HEADCROP:10"],
      extra=""
   threads: 8
   wrapper:
      "0.27.1/bio/trimmomatic/se"
