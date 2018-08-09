#!/bin/bash -l

# load uppmax modules
module load slurm/16.05.11
module load samtools/1.3.1
module load blast/2.2.25
module load oracle-jdk/1.7.0_79
module load perl/5.24.0
module load perl/bioperl/1.6.1
module load gnuplot/5.0.2
module load fastqc/0.11.5
module load cutadapt/1.8.0
module load bwa/0.7.12
module load trimmomatic/0.35
module load parallel/20161022
module load pindel/0.2.5a8
module load trimmomatic/0.35

export PYTHONPATH="/home/patsm159/workspace/CAPS:/home/patsm159/anaconda3/lib/python3.6/site-packages"
 snakemake  -j 128 --cluster-config /home/patsm159/workspace/CAPS/cluster.json --cluster "sbatch -A wp1 -p core -n {cluster.n} -t 24:00:00 --output=logs/slurm/slurm_%j.out" --directory /projects/wp1/nobackup/workspace/test_snakemake/ -s /projects/wp1/nobackup/workspace/test_snakemake/Snakefile --verbose --wrapper-prefix git+file:///home/patsm159/workspace/merged-snakemakewrappers
