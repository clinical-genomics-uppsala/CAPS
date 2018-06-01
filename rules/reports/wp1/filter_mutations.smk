# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.data.report.wp1 import generate_filtered_mutations

def get_sample_value(field, sample, default = False, converter = None):
    if field in samples and sample in samples[field]:
        value = samples[field][sample]
        if converter is None:
            return value
        else:
            return converter(value)
    else:
        return default

def converter(value):
    if value == "yes":
        return True
    elif value == "no":
        return False
    else:
        raise Exception("Unexpected value: " + str(value) + ", should have been yes or no.")

rule filter_mutations:
    input:
        snpmania = 'jsnpmania/{sample}.variations',
        annovar = ['pindel_annovar_output/{sample}.pindel.singleSample.annovarOutput', 'annovar_output/{sample}.singleSample.annovarOutput']
    output:
        "reports/{sample}.filteredMutations.tsv"
    params:
        amplicon_mapped = lambda wildcards: get_sample_value('amplicon_mapped', wildcards.sample, converter=converter),
        blacklist = lambda wildcards: config["files"][samples['panel_type'][wildcards.sample]]['blacklist'],
        chr_to_nc = config["chr_to_nc_file"],
        hotspot = lambda wildcards: samples["hotspot"][wildcards.sample],
        max_1000genome = lambda wildcards: config['filter_settings'][samples['panel_type'][wildcards.sample]]['max_1000genome_ratio'],
        min_read_depth = lambda wildcards: config['filter_settings'][samples['panel_type'][wildcards.sample]]['min_read_depth'],
        min_vaf = lambda wildcards: config['filter_settings'][samples['panel_type'][wildcards.sample]]['min_allele_ratio'],
        multibp = config["files"]["multibp_variants_mappers"],
        read_depth_classes = lambda wildcards: [(300 ,"ok","yes"), (30 ,"low","yes"),(0 ,"lowk","not analyzable")],
        transcript = config["files"]["transcripts_mappers"]
    log:
        "logs/reports/wp1/{sample}.filter_mutations.log"
    run:
        generate_filtered_mutations(
            wildcards.sample,
            output[0],
            params.hotspot,
            input.snpmania,
            input.annovar,
            params.multibp,
            params.chr_to_nc,
            params.min_read_depth,
            params.min_vaf,
            params.read_depth_classes,
            params.max_1000genome,
            params.amplicon_mapped,
            params.blacklist,
            params.transcript)
