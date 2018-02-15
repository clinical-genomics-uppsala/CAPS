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
        temp("pindel/{sample}.indels_D.temp.vcf")
    params:
        refname=config['reference_genome_name'],
        refdate=['reference_genome_date'],
        extra="--min_size 3"
    log:
        "logs/pindel/pindel2vcf.D.log"
    wrapper:
        "0.19.3/bio/pindel/pindel2vcf"

rule check_vcf_D_for_samplename:
    input:
        "pindel/{sample}.indels_D.temp.vcf"
    output:
        "pindel/{sample}.indels_D.vcf"
    shell:
        "sed 's/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO$/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t{wildcards.sample}/' {input} > {output}"

#ToDo Remember update wrapper after pull-request has been approved
rule pindel_to_vcf:
    input:
        ref=config['reference_genome'],
        pindel="pindel/{sample}.indels_D"
    output:
        "pindel/{sample}.vcf"
    params:
        refname=config['reference_genome_name'],
        refdate=['reference_genome_date'],
        extra="--min_size 3",
        prefix="pindel/{sample}.indels"
    log:
        "logs/pindel/pindel2vcf.SI.log"
    wrapper:
        "https://bitbucket.org/Smeds/snakemake-wrappers/raw/a2341cdc21018dd43a2078aa542119fcdb89bb7c/bio/pindel/pindel2vcf/wrapper.py"
