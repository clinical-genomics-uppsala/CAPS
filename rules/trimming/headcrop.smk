# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq

#rule headcrop_read:
#   input:
#      lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
#   output:
#      temp("trimmed/{sample}.{unit}.{read}.headcrop.fastq.gz")
#   params:
#      trimmer=["HEADCROP:10"],
#      extra=""
#   threads: 8
#   wrapper:
#      "0.27.1/bio/trimmomatic/se"

rule headcrop_read_part:
   input:
        "trimmed/{sample}.{unit}.{part}.{read}.fastq.gz"
   output:
      temp("trimmed/{sample}.{unit}.{part}.{read}.headcrop.fastq.gz")
   log:
      "logs/trimming/{sample}.{unit}.{part}.{read}.headcrop.trimmomatic.txt"
   params:
      trimmer=["HEADCROP:10"],
      extra=""
   threads: 8
   wrapper:
      "0.27.1/bio/trimmomatic/se"
