import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index


qc_files = [("fastqc", ".html"),("fastqc" ,".zip")]
def generate_file_output_qc():
    return [f[0] + "/" + str(row.Index[0]) + "." + str(row.Index[1]) + f[1]
        for row in units.itertuples()
            for f in qc_files]

def generate_bam():
    return [os.path.join("calls/", str(row.Index) + ".lofreq.merged.vcf") for row in samples.itertuples()]
rule all:
    input:
        generate_bam()

include: "workflows/accel_amplicon_hs.snakemake"
