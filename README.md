# CAPS : Clinical Analysis Pipeline - Snakemae

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/clinical-genomics-uppsala/CAPS.svg?branch=master)](https://travis-ci.org/clinical-genomics-uppsala/CAPS)

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description, covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules` and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

## Authors

* Patrik Smeds (@smeds)

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/caps/releases).
If you intend to modify and further develop this workflow, fork this reposity. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.


snakemake  -j 128 --cluster-config cluster.json --cluster "sbatch -A wp1 -p core -n {cluster.n} -t 24:00:00 --output=logs/slurm/slurm_%j.out" --directory /projects/wp1/nobackup/workspace/ -s Snakefile --verbose --wrapper-prefix git+file:///home/patsm159/workspace/merged-snakemakewrappers --latency-wait 30 --restart-times 10

SINGULARITYENV_PREPEND_PATH="$PATH:/bianca/bin" singularity exec -B /var/run/munge/:/run/munge -B /usr/bin/:/bianca/bin/ -B /usr/lib64,/etc/slurm,/etc/passwd ../singularity/caps.simg snakemake -j 128 --cluster-config /proj/sens2017561/nobackup/wharf/patriksm/patriksm-sens2017561/CAPS/cluster.json --cluster "sbatch -A sens2017561 -n {cluster.n} -t 24:00:00 --output=logs/slurm/slurm_%j.out" -s /proj/sens2017561/nobackup/wharf/patriksm/patriksm-sens2017561/CAPS/Snakefile --verbose --wrapper-prefix git+file:///castor/project/proj_nobackup/wharf/patriksm/patriksm-sens2017561/snakemake-wrappers-source --verbose --js jobscript_mod.sh
 1000  cd /proj/sens2017561/nobackup/wharf/patriksm/patriksm-sens2017561/

## Testing

Tests cases are in the subfolder `.test`. They should be executed via continuous integration with Travis CI.

## Local Testing

docker run -v /Users/patsm159/auctornotitia_workspace/caps-snakemake:/snakemake-workflows -it snakemake bash

'''
pip install pyvcf
pip install biopython
apt update
apt install -y ttf-dejavu
apt install -y build-essential
apt install -y libpod-plainer-perl
'''
