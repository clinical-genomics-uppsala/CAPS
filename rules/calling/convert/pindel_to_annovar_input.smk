# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#from scripts.lib.common.data.parser.pindel import convert_to_annovar_input

rule pindel_to_annovar_input:
    input:
        pindel_deletions="pindel/{sample}.indels_D",
        pindel_insertions="pindel/{sample}.indels_SI",
        vcf="pindel/{sample}.merged.vcf"
    output:
        "PindelAnnovar/{sample}.pindel.filtered.annovarInput"
    #params:
    #    min_read_depth = lambda wildcards: config["pindel_flags"][samples["panel_type"][wildcards.sample]]['min_read_depth'],
    #    min_allele_ratio = lambda wildcards: config["pindel_flags"][samples["panel_type"][wildcards.sample]]['min_allele_ratio'],
    #    read_method = lambda wildcards: config["pindel_flags"][samples["panel_type"][wildcards.sample]]['read_method']
    #log:
#        "logs/jsnpmania/{sample}.jsnpmania_to_annovar.log"
    run:
        print("TESt")
