# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#pindel_types = ["D", "BP", "INV", "TD", "LI", "SI", "RP"]
pindel_types = ["D", "SI"]

from scripts.lib.common.utils import get_bam_file

rule create_pindel_config:
    input:
        lambda wildcards: "mapped/" + get_bam_file(wildcards, samples, True)
    output:
        "pindel/{sample}.config.txt"
    shell:
        "echo -e \"{input}\t300\t{wildcards.sample}\" > {output}"

rule pindel:
    input:
        ref=config['reference_genome'],
        #samples=lambda wildcards: "mapped/" + get_bam_file(wildcards, samples, True),
        samples="mapped/{sample}.sorted.bam",
        config="pindel/{sample}.config.txt",
        index=lambda wildcards: "mapped/" + get_bam_file(wildcards, samples,True) + ".bai"
    output:
        expand("pindel/{{sample}}.indels_{type}", type=pindel_types)
    params:
        prefix= "pindel/{sample}.indels",# lambda wildcards: "pindel/" + wildcards.sample + ".indels",
        extra= lambda wildcards: " -J " + config['pindel_exclude_regions'] + \
            " -x " + str(config['pindel_flags']['max_range_index']) + \
            " -B " + str(config['pindel_flags']['balance_cutoff'])
    log:
        "logs/pindel/{sample}.log"
    threads: 5
    wrapper:
        "0.24.0/bio/pindel/call"

rule pindel_to_vcf:
    input:
        ref=config['reference_genome'],
        pindel=["pindel/{sample}.indels_D", "pindel/{sample}.indels_SI"] # Waiting for pull-request to be approved
    output:
        "pindel/{sample}.vcf"
    params:
        refname=config['reference_genome_name'],
        refdate=config['reference_genome_date'],
        extra="--min_size 3"
    log:
        "logs/pindel/{sample}.pindel2vcf.log"
    wrapper:
       "0.26.1/bio/pindel/pindel2vcf"
