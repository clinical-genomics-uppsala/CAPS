# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

#file_endings={"Swift": [".R1.trimmomatic_cutadapt.fastq.gz", ".R2.trimmomatic_cutadapt.fastq.gz", "_001.trimmomatic_cutadapt.qc.txt"],
#              "HaloPlex": [".cutadapt.fastq.gz", "_R2_001.cutadapt.fastq.gz", "_001.cutadapt.qc.txt"]}

#file_endings = [".variations", ".deletions", ".insertions"]
#file_endings = [".cutadapt.fastq.gz", "_R2_001.cutadapt.fastq.gz", "_001.cutadapt.qc.txt"]
#pindel_file_endings = [".indels_D",".indels_SI"]
pindel_file_endings = [".vcf"]
def generate_file_output_pindel():
    return [os.path.join("pindel", str(row.Index) + ending) for row in samples.itertuples() for ending in pindel_file_endings]

qc_files = [("fastqc", ".html"),("fastqc" ,".zip")]
def generate_file_output_qc():
    return [f[0] + "/" + str(row.Index[0]) + "-" + str(row.Index[1]) + f[1]
        for row in units.itertuples()
            for f in qc_files]

def generate_file_output_jsnpmania_input():
    return [os.path.join("jsnpmania", str(row.Index) + ending) for row in samples.itertuples() for ending in ['.variations','.insertions','.deletions']]

def generate_file_output_annovar_input():
    return [os.path.join("annovar", str(row.Index) + ".annovarInput") for row in samples.itertuples()]

def generate_pindel_file_output_annovar_input():
    return [os.path.join("pindel_annovar", str(row.Index) + ".pindel.filtered.annovarInput") for row in samples.itertuples()]

def generate_filtered_mutations():
    return [os.path.join("reports", str(row.Index) + ".filteredMutations.tsv") for row in samples.itertuples()]
rule all:
    input:
        generate_file_output_qc() + generate_pindel_file_output_annovar_input()

include: "workflows/wp1.snakemake"
