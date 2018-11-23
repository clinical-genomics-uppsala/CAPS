# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_fastq

def _cgu_get_num_splits(config):
    return int(config.get("num_fastq_split",1))

rule count_lines_in_fastq:
    input:
      lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
    output:
      temp("trimmed/{sample}.{unit}.{read}.var")
    run:
      import subprocess, os
      lines = int(float(subprocess.run("gunzip -c " + input[0] + " |  wc -l | awk '{print $1/4}'", stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8').rstrip("\n")))
      storage.store(wildcards.sample + "." + wildcards.unit + "." + wildcards.read + ".var",str(lines))
      shell("echo 'reads: '" + str(lines) + "'' > "  + output[0])

rule split_fastq_file:
    input:
      lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2'),
      "trimmed/{sample}.{unit}.{read}.var"
    output:
      temp(['trimmed/{sample}.{unit}.%04d.{read}.fastq' % num for num in range(0,_cgu_get_num_splits(config))])
    params:
      output_prefix=lambda wildcards: "trimmed/" + wildcards.sample + "." + wildcards.unit + ".",
      output_suffix=lambda wildcards: "." + wildcards.read + ".fastq"
    run:
      import math
      num_reads = int(storage.fetch(wildcards.sample + "." + wildcards.unit + "." + wildcards.read + ".var"))
      num_split = _cgu_get_num_splits(config)
      lines_per_file = 4*math.ceil(num_reads / num_split)
      shell('gunzip -c {input[0]} | awk \'BEGIN{{ file = 0; filename = sprintf("{params.output_prefix}%.04d{params.output_suffix}", file) }}{{ print > filename}} NR % {lines_per_file} == 0 {{ close(filename); file = file + 1; filename = sprintf("{params.output_prefix}%.04d{params.output_suffix}",file)}}\'')
      num_files_generated = 4*math.floor(num_reads / lines_per_file)
      while num_files_generated < num_split:
        shell("touch {params.output_prefix}%04d{params.output_suffix}" % num_split)
        num_split -= 1

rule compress_split_fastq:
    input:
        "trimmed/{sample}.{unit}.{part}.{read}.fastq"
    output:
        "trimmed/{sample}.{unit}.{part}.{read}.fastq.gz"
    run:
        shell("gzip -c {input} > {output}")
