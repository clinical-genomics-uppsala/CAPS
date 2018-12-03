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

storage = PersistentDict("caps_mystorage")

_cutadapt_input = ["trimmed/{sample}.{unit}.{part}.R1.fastq", "trimmed/{sample}.{unit}.{part}.R2.fastq"]
try:
    _cutadapt_input = cutadapt_input
except:
    pass

_cutadapt_fastq1_output = "trimmed/{sample}.{unit}.{part}.R1.cutadapt.fastq.gz"
try:
    _cutadapt_fastq1_output = cutadapt_fastq1_output
except:
    pass

_cutadapt_fastq2_output = "trimmed/{sample}.{unit}.{part}.R2.cutadapt.fastq.gz"
try:
    _cutadapt_fastq2_output = cutadapt_fastq2_output
except:
    pass

rule cutadapt:
  input:
      _cutadapt_input
  output:
      fastq1=_cutadapt_fastq1_output,
      fastq2=_cutadapt_fastq2_output,
      qc = "logs/trimmed/{sample}.{unit}.{part}.cutadapt.qc.txt"
  params:
      " --minimum-length 1",
      lambda wildcards: " -a " + config["cutadapt"][samples['adapters'][sample_id(wildcards)]]["first_pair_adapter"],
      lambda wildcards: " -A " + reverse_complement(config["cutadapt"][samples['adapters'][sample_id(wildcards)]]["second_pair_adapter"])
  wrapper:
      "0.17.4/bio/cutadapt/pe"
