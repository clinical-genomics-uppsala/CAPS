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

## Testing

Tests cases are in the subfolder `.test`. They should be executed via continuous integration with Travis CI.

## Local Testing

docker run -v /Users/patsm159/auctornotitia_workspace/caps-snakemake:/snakemake-workflows -it snakemake bash

'''
apt update
apt install ttf-dejavu
apt-get install build-essential
apt install libpod-plainer-perl
'''
