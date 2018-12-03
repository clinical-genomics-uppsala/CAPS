# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.data.report.wp1 import generate_filtered_mutations, filtered_mutations_to_vcf


def field_converter(value):
    if value == "yes":
        return True
    elif value == "no":
        return False
    else:
        raise Exception("Unhandled field value: " + value)

def get_sample_value(field, sample, converter=None, default=False):
    if field in samples and sample in samples[field]:
        return converter(samples[field][sample])
    else:
        return default

rule filter_mutations:
    input:
        snpmania = 'jsnpmania/{sample}.variations',
        annovar = ['annovar_output/{sample}.pindelFiltered.singleSample.annovarOutput', 'annovar_output/{sample}.jsnpmania.singleSample.annovarOutput']
    output:
        "reports/{sample}.filteredMutations.tsv"
    params:
        amplicon_mapped = lambda wildcards: get_sample_value('amplicon_mapped', wildcards.sample, converter=field_converter),
        blacklist = lambda wildcards: config["files"][samples['panel_type'][wildcards.sample]]['blacklist'],
        chr_to_nc = config["chr_to_nc_file"],
        hotspot = lambda wildcards: samples["hotspot"][wildcards.sample],
        max_1000genome = lambda wildcards: config['filter_settings'][samples['sample_source'][wildcards.sample]]['max_1000genome_ratio'],
        min_read_depth = lambda wildcards: config['filter_settings'][samples['sample_source'][wildcards.sample]]['min_read_depth'],
        min_vaf = lambda wildcards: config['filter_settings'][samples['sample_source'][wildcards.sample]]['min_allele_ratio'],
        multibp = config["files"]["multibp_variants_mappers"],
        read_depth_classes = lambda wildcards: [(300 ,"ok","yes"), (30 ,"low","yes"),(0 ,"low","not analyzable")],
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

rule filter_mutations_to_vcf:
  input:
    "reports/{sample}.filteredMutations.tsv"
  output:
    "reports/{sample}.vcf"
  params:
    ampregion_file = lambda wildcards: samples["snpseq_file"][wildcards.sample],
    chr_to_nc = config["chr_to_nc_file"]
  log:
      "logs/reports/wp1/{sample}.filter_mutations_to_vcf.log"
  run:
    filtered_mutations_to_vcf(
      input[0],
      output[0],
      params.ampregion_file,
      params.chr_to_nc)

rule filter_mutations_to_vcf_vaf0_05:
  input:
    "reports/{sample}.filteredMutations.tsv"
  output:
    "reports/{sample}.vaf0.05.vcf"
  params:
    ampregion_file = lambda wildcards: samples["snpseq_file"][wildcards.sample],
    chr_to_nc = config["chr_to_nc_file"]
  log:
      "logs/reports/wp1/{sample}.filter_mutations_to_vcf_vaf0_05.log"
  run:
    filtered_mutations_to_vcf(
      input[0],
      output[0],
      params.ampregion_file,
      params.chr_to_nc,
      min_vaf=0.05)
