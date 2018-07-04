# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""

Expects the global variable config of at least the following structure
...............................................................................
---
illuminaclip_file: /path/to/illumina.fa
cutadapt:
  adapters_1:
    first_pair_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC,
    second_pair_adapter: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
...............................................................................

Expects a samples config
With the following dictionary structure (panel, sample_number and lane are
required)
samples = {
    'cutadapt': {sample_id: value},
    'fq1': {samplei_id: value},
    'fq2' {sample_id: value}
}

Example of a sample.tsv that can be imported using pandas (columns need to be
tab separated)
...............................................................................
sample     adapters     fq1                fq2
sample1    adapters_1   sample1.R1.fastq   sample1.R2.fastq
sample2    adapters_1   sample2.R1.fastq   sample2.R2.fastq
...............................................................................
"""
from scripts.lib.common.utils import get_fastq, sample_id, reverse_complement
from pytools.persistent_dict import PersistentDict

storage = PersistentDict("mystorage")

rule extract_fastq_files:
   input:
      lambda wildcards: get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
   output:
      temp("trimmed/.temp/{sample}.{unit}.{read}.fastq")
   run:
      if input[0].endswith("gz"):
         shell("zcat {input} > {output}")
      else:
         shell("cat {input} > {output}")


rule count_lines_in_fastq:
    input:
      "trimmed/.temp/{sample}.{unit}.{read}.fastq"
    output:
      temp("trimmed/.temp/{sample}.{unit}.{read}.var")
    run:
      import subprocess, os
      lines = int(float(subprocess.run("wc -l " + str(input[0]) + " |  awk '{print($1/4)}'", stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8').rstrip("\n")))
      storage.store(wildcards.sample + "." + wildcards.unit + "." + wildcards.read + ".var",str(lines))
      shell("echo 'reads: '" + str(lines) + "'' > "  + output[0])

rule split_fastq_file:
   input:
      "trimmed/.temp/{sample}.{unit}.{read}.fastq",
      "trimmed/.temp/{sample}.{unit}.{read}.var"
   output:
      temp(['trimmed/.temp/{sample}.{unit}.%02d.{read}.fastq' % num for num in range(0,int(config.get("num_fastq_split",1)))])
   params:
      output_prefix=lambda wildcards: "trimmed/.temp/" + wildcards.sample + "." + wildcards.unit + ".",
      output_suffix=lambda wildcards: "." + wildcards.read + ".fastq"
   run:
     num_reads = int(storage.fetch(wildcards.sample + "." + wildcards.unit + "." + wildcards.read + ".var"))
     num_split = int(config.get("num_fastq_split",1))
     import math
     lines_per_file = 4*math.ceil(num_reads / num_split)
     shell("split -d -l {lines_per_file} {input[0]} {params.output_prefix} --additional-suffix={params.output_suffix}")
     if num_reads < lines_per_file*num_split:
        number_of_generated_files = num_split - math.ceil(num_reads/num_split)
        if num_split > number_of_generated_files:
            for i in range(number_of_generated_files,num_split):
                shell("touch {params.output_prefix}%02d{params.output_suffix}" % i)


rule cutadapt:
   input:
       [lambda wildcards: get_fastq(wildcards,units,'fq1'),
        lambda wildcards: get_fastq(wildcards,units,'fq2')]
   output:
       fastq1="trimmed/{sample}.{unit}.R1.cutadapt.fastq.gz",
       fastq2="trimmed/{sample}.{unit}.R2.cutadapt.fastq.gz",
       qc = "logs/trimmed/{sample}.{unit}.cutadapt.qc.txt"
   params:
       " --minimum-length 1",
       lambda wildcards: " -a " + config["cutadapt"][samples['adapters'][sample_id(wildcards)]]["first_pair_adapter"],
       lambda wildcards: " -A " + reverse_complement(config["cutadapt"][samples['adapters'][sample_id(wildcards)]]["first_pair_adapter"])
   wrapper:
       "0.17.4/bio/cutadapt/pe"

rule cutadapt_split:
  input:
      ["trimmed/.temp/{sample}.{unit}.{part}.R1.fastq",
       "trimmed/.temp/{sample}.{unit}.{part}.R2.fastq"]
  output:
      fastq1="trimmed/{sample}.{unit}.{part}.R1.cutadapt.fastq.gz",
      fastq2="trimmed/{sample}.{unit}.{part}.R2.cutadapt.fastq.gz",
      qc = "logs/trimmed/{sample}.{unit}.{part}.cutadapt.qc.txt"
  params:
      " --minimum-length 1",
      lambda wildcards: " -a " + config["cutadapt"][samples['adapters'][sample_id(wildcards)]]["first_pair_adapter"],
      lambda wildcards: " -A " + reverse_complement(config["cutadapt"][samples['adapters'][sample_id(wildcards)]]["first_pair_adapter"])
  wrapper:
      "0.17.4/bio/cutadapt/pe"

#rule merge_cutadapt
#  input:
#    lambda wildcards: ["trimmed/" + wildcards.sample + "." + wildcards.unit + "." + p + "." + wildcards.read
#  output:
#    fastq1="trimmed/.temp/{sample}_{unit}_R1/{sample}.{unit}.{read}.cutadapt.fastq.gz"
