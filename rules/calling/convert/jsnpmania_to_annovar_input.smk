# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.data.parser.jsnpmania import convert_jsnpmania_to_annvovar_output


_jsnpmania_to_annovar_input_variations_input = "jsnpmania/{sample}.variations"
try:
    _jsnpmania_to_annovar_input_variations_input = jsnpmania_to_annovar_input_variations_input
except:
    pass

_jsnpmania_to_annovar_input_deletions_input = "jsnpmania/{sample}.deletions"
try:
    _jsnpmania_to_annovar_input_deletions_input = jsnpmania_to_annovar_input_deletions_input
except:
    pass

_jsnpmania_to_annovar_input_inseriton_input = "jsnpmania/{sample}.insertions"
try:
    _jsnpmania_to_annovar_input_inseriton_input = jsnpmania_to_annovar_input_inseriton_input
except:
    pass

_jsnpmania_to_annovar_input_output = "annovar/{sample}.jsnpmania.annovarInput"
try:
    _jsnpmania_to_annovar_input_output = jsnpmania_to_annovar_input_output
except:
    pass

rule jsnpmania_to_annovar_input:
    input:
        jsnpmania_variants=_jsnpmania_to_annovar_input_variations_input,
        jsnpmania_insertions=_jsnpmania_to_annovar_input_inseriton_input,
        jsnpmania_deletions=_jsnpmania_to_annovar_input_deletions_input
    output:
        _jsnpmania_to_annovar_input_output
    params:
        nc_to_chr_file = config['nc_to_chr_file'],
        min_allele_ratio = lambda wildcards: config["filter_settings"][samples["sample_source"][wildcards.sample]]['min_allele_ratio'],
        min_read_depth = lambda wildcards: config["filter_settings"][samples["sample_source"][wildcards.sample]]['min_read_depth_snpmania'],
        amplicon_min_depth = lambda wildcards: config["filter_settings"][samples["sample_source"][wildcards.sample]]['amplicon_min_depth']
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
