# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.data.parser.jsnpmania import convert_jsnpmania_to_annvovar_output

rule jsnpmania_to_annovar_input:
    input:
        jsnpmania_variants="JSNPmania/{sample}.variations",
        jsnpmania_insertions="JSNPmania/{sample}.insertions",
        jsnpmania_deletions="JSNPmania/{sample}.deletions"
    output:
        "Annovar/{sample}.annovarInput"
    params:
        nc_to_chr_file = config['nc_to_chr_file'],
        min_allele_ratio = lambda wildcards: config["filter_settings"][samples["panel_type"][wildcards.sample]]['min_allele_ratio'],
        min_read_depth = lambda wildcards: config["filter_settings"][samples["panel_type"][wildcards.sample]]['min_read_depth'],
        amplicon_min_depth = lambda wildcards: config["filter_settings"][samples["panel_type"][wildcards.sample]]['amplicon_min_depth']
    log:
        "logs/jsnpmania/{sample}.jsnpmania_t_annovar.log"
    run:
        convert_jsnpmania_to_annvovar_output(
            wildcards.sample,
            output[0],
            input.jsnpmania_variants,
            input.jsnpmania_insertions,
            input.jsnpmania_deletions,
            params.nc_to_chr_file,
            params.min_allele_ratio,
            params.min_read_depth,
            params.amplicon_min_depth)
