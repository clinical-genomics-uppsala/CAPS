# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#pindel_types = ["D", "BP", "INV", "TD", "LI", "SI", "RP"]
pindel_types = ["D", "SI"]

from scripts.lib.common.utils import get_bam_file

_pindel_input = "mapped/{sample}.sorted.bam"
try:
    _pindel_input = pindel_input
except:
    pass

_pindel_output = "pindel/{{sample}}.indels_{type}"
try:
    _pindel_output = pindel_output
except:
    pass

_pindel_vcf_output = "pindel/{sample}.vcf"
try:
    _pindel_vcf_output = pindel_vcf_output
except:
    pass

rule create_pindel_config:
    input:
        _pindel_input
    output:
        temp("pindel/{sample}.config.txt")
    shell:
        "echo -e \"{input}\t300\t{wildcards.sample}\" > {output}"

rule pindel:
    input:
        ref=config['reference_genome'],
        #samples=lambda wildcards: "mapped/" + get_bam_file(wildcards, samples, True),
        samples=_pindel_input,
        config="pindel/{sample}.config.txt",
        index=_pindel_input + ".bai"
    output:
        expand(_pindel_output, type=pindel_types)
    params:
        prefix= "pindel/{sample}.indels",# lambda wildcards: "pindel/" + wildcards.sample + ".indels",
        extra= lambda wildcards: " -J " + config['pindel_exclude_regions'] + \
            " -x " + str(config['pindel_flags']['max_range_index']) + \
            " -B " + str(config['pindel_flags']['balance_cutoff'])
    log:
        "logs/pindel/{sample}.log"
    threads: 5
    wrapper:
        "0.27.1//bio/pindel/call"

rule pindel_to_vcf:
    input:
        ref=config['reference_genome'],
        pindel=expand(_pindel_output, type=pindel_types) # Waiting for pull-request to be approved
    output:
        _pindel_vcf_output
    params:
        refname=config['reference_genome_name'],
        refdate=config['reference_genome_date'],
        extra="--min_size 3"
    log:
        "logs/pindel/{sample}.pindel2vcf.log"
    wrapper:
       "0.27.1//bio/pindel/pindel2vcf"
