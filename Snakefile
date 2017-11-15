# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

#file_endings={"Swift": [".R1.trimmomatic_cutadapt.fastq.gz", ".R2.trimmomatic_cutadapt.fastq.gz", "_001.trimmomatic_cutadapt.qc.txt"],
#              "HaloPlex": [".cutadapt.fastq.gz", "_R2_001.cutadapt.fastq.gz", "_001.cutadapt.qc.txt"]}

file_endings = [".sorted.bam",".sorted.bam.bai"]
#file_endings = [".cutadapt.fastq.gz", "_R2_001.cutadapt.fastq.gz", "_001.cutadapt.qc.txt"]
def generate_file_output():
    return [os.path.join("mapped", str(row.Index[0]) + "-" + str(row.Index[1]) + ending) for row in units.itertuples() for ending in file_endings]


rule all:
    input:
        generate_file_output()

print(generate_file_output())
include: "workflows/wp1.snakemake"
