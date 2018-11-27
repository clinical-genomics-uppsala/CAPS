from scripts.lib.common.utils import get_fastq

_extract_umis_input = lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
try:
    _extract_umis_input = extract_umis_input
except:
    pass

_extract_umis_output = temp("umi/{sample}.{unit}.{part}.{read}.UMIs.fastq")
try:
    _extract_umis_output = extract_umis_output
except:
    pass

rule extract_UMIs:
    input:
       _extract_umis_input
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
      _extract_umis_output
   params:
      trimmer=["CROP:10"],
      extra=""
   threads: 8
   wrapper:
      "0.27.1/bio/trimmomatic/se"
