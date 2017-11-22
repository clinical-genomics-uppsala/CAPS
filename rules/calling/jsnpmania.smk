# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.common.utils import get_bam_file

rule run_jsnpmania:
    input:
        lambda wildcards: "mapped/" + get_bam_file(wildcards, samples)
    output:
        variations="JSNPmania/{sample}.variations",
        insertions="JSNPmania/{sample}.insertions",
        deletions="JSNPmania/{sample}.deletions"
    params:
        path_jsnpmania = config.get('path_jsnpmania',"jsnpmania.jar"),
        path_jsnpmania_header = config['path_jsnpmania_header'],
        ref_file=lambda wildcards: samples['snpseq_file'][wildcards.sample],
        flags= lambda wildcards: config["jsnpmania_flags"]
    log:
        "logs/jsnpmania/{sample}.jsnpmania.log"
    wrapper:
        "https://raw.githubusercontent.com/clinical-genomics-uppsala/snakemake-wrappers/add-java-dep/bio/jsnpmania/wrapper.py"
