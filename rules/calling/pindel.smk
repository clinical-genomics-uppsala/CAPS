# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#pindel_types = ["D", "BP", "INV", "TD", "LI", "SI", "RP"]
pindel_types = ["D"]

from scripts.common.utils import get_bam_file

rule create_pindel_config:
    input:
        lambda wildcards: "mapped/" + get_bam_file(wildcards, samples)
    output:
        "pindel/{sample}.config.txt"
    shell:
        "echo -e \"{input}\t300\t{wildcards.sample}\" > {output}"

rule pindel:
    input:
        ref=config['reference_genome'],
        samples=lambda wildcards: "mapped/" + get_bam_file(wildcards, samples),
        config="pindel/{sample}.config.txt",
        index=lambda wildcards: "mapped/" + get_bam_file(wildcards, samples) + ".bai"
    output:
        expand("pindel/{{sample}}_{type}", type=pindel_types)
    params:
        prefix="pindel/{sample}",
        extra= lambda wildcards: " -J " + config['pindel_exclude_regions'] + \
            " " + config["pindel_flags"][samples["panel_type"][wildcards.sample]]
    log:
        "logs/pindel/{sample}.log"
    threads: 4
    wrapper:
        "0.19.3/bio/pindel/call"
