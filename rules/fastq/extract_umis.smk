# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq

rule extract_UMIs_head_part:
   input:
      "trimmed/{sample}.{unit}.{part}.{read}.fastq.gz"
   output:
      temp("umi/{sample}.{unit}.{part}.{read}.UMIs.fastq")
   log:
      "logs/umi/{sample}.{unit}.{part}.{read}.UMI.head.txt"
   params:
      trimmer=["CROP:10"],
      extra=""
   threads: 8
   wrapper:
      "0.27.1/bio/trimmomatic/se"
