# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.data.parser.pindel import convert_to_annovar_input

_pindel_to_annovar_d_input = "pindel/{sample}.indels_D"
try:
    _pindel_to_annovar_d_input = pindel_to_annovar_d_input
except:
    pass

_pindel_to_annovar_i_input = "pindel/{sample}.indels_SI"
try:
    _pindel_to_annovar_i_input = pindel_to_annovar_i_input
except:
    pass

_pindel_to_annovar_i_output = "annovar/{sample}.pindelFiltered.annovarInput"
try:
    _pindel_to_annovar_i_output = pindel_to_annovar_i_output
except:
    pass

rule pindel_to_annovar_input:
    input:
        pindel_deletions=_pindel_to_annovar_d_input,
        pindel_insertions=_pindel_to_annovar_i_input,
        vcf="pindel/{sample}.vcf"
    output:
        _pindel_to_annovar_i_output
    params:
        min_read_depth = lambda wildcards: config["annovar_pindel_flags"][samples["sample_source"][wildcards.sample]]['min_read_depth'],
        min_allele_ratio = lambda wildcards: config["annovar_pindel_flags"][samples["sample_source"][wildcards.sample]]['min_allele_ratio'],
        read_method = lambda wildcards: config["annovar_pindel_flags"][samples["sample_source"][wildcards.sample]]['read_method']
    log:
        "logs/pindel/{sample}.pindel_to_annovar.log"
    run:
        convert_to_annovar_input(
            wildcards.sample,
            output[0],
            input.vcf,
            input.pindel_deletions,
            input.pindel_insertions,
            params.min_read_depth,
            params.min_allele_ratio,
            params.read_method)
