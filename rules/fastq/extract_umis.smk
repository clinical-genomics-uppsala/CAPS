# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq

_extract_umis_head_input = "trimmed/{sample}.{unit}.{part}.{read}.fastq.gz"
try:
    _extract_umis_head_input = extract_umis_head_input
except:
    pass

_extract_umis_head_output = temp("umi/{sample}.{unit}.{part}.{read}.UMIs.fastq")
try:
    _extract_umis_head_output = extract_umis_head_output
except:
    pass

rule extract_UMIs_head:
   input:
      _extract_umis_head_input
   output:
      _extract_umis_head_output
   log:
      "logs/umi/{sample}.{unit}.{part}.{read}.UMI.head.txt"
   params:
      trimmer=["CROP:10"],
      extra=""
   threads: 8
   wrapper:
      "0.27.1/bio/trimmomatic/se"
