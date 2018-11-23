from scripts.lib.common.utils import get_fastq

rule extract_UMIs:
    input:
       lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
    output:
       temp("umi/{sample}.{unit}.{read}.UMIs.fastq")
    params:
      trimmer=["CROP:10"],
       extra=""
    threads: 8
    wrapper:
       "0.27.1/bio/trimmomatic/se"

rule extract_UMIs_part:
   input:
      "trimmed/{sample}.{unit}.{part}.{read}.fastq.gz"
   output:
      temp("umi/{sample}.{unit}.{part}.{read}.UMIs.fastq")
   params:
      trimmer=["CROP:10"],
      extra=""
   threads: 8
   wrapper:
      "0.27.1/bio/trimmomatic/se"
