# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#pindel_types = ["D", "BP", "INV", "TD", "LI", "SI", "RP"]
pindel_types = ["D", "SI"]

from scripts.lib.common.utils import get_bam_file

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
        expand("pindel/{{sample}}.indels_{type}", type=pindel_types)
    params:
        prefix= "pindel/{sample}.indels",# lambda wildcards: "pindel/" + wildcards.sample + ".indels",
        extra= lambda wildcards: " -J " + config['pindel_exclude_regions'] + \
            " -m " + str(config["pindel_flags"][samples["panel_type"][wildcards.sample]]['min_read_depth']) +  \
            " -v " + str(config["pindel_flags"][samples["panel_type"][wildcards.sample]]['min_allele_ratio']) +  \
            " --report_long_insertions False " +  \
            " --report_close_mapped_reads False" +  \
            " --report_breakpoints False" +  \
            " --report_inversions False" +  \
            " --report_duplications False"
    log:
        "logs/pindel/{sample}.log"
    threads: 4
    wrapper:
        "0.19.3/bio/pindel/call"

rule pindel_to_vcf_D:
    input:
        ref=config['reference_genome'],
        pindel="pindel/{sample}.indels_D"
    output:
        "pindel/{sample}.indels_D.vcf"
    params:
        refname=config['reference_genome_name'],
        refdate=['reference_genome_date'],
        extra="--min_size 3"
    log:
        "logs/pindel/pindel2vcf.D.log"
    wrapper:
        "0.19.3/bio/pindel/pindel2vcf"

rule pindel_to_vcf_SI:
    input:
        ref=config['reference_genome'],
        pindel="pindel/{sample}.indels_SI"
    output:
        "pindel/{sample}.indels_SI.vcf"
    params:
        refname=config['reference_genome_name'],
        refdate=['reference_genome_date'],
        extra="--min_size 3"
    log:
        "logs/pindel/pindel2vcf.SI.log"
    wrapper:
        "0.19.3/bio/pindel/pindel2vcf"

rule merge_pindel_vcf_files:
    input:
        calls=["pindel/{sample}.indels_D.vcf","pindel/{sample}.indels_SI.vcf"]
    output:
        "pindel/{sample}.merged.vcf"
    log:
        "logs/pindel/bcftools.merge.log"
    wrapper:
        "0.19.3/bio/bcftools/concat"
