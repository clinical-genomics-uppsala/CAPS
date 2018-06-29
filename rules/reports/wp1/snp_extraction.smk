# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.data.parser.jsnpmania import extract_snp_file, position_filter

rule extract_snp:
    input:
        "jsnpmania/{sample}.variations"
    output:
        "reports/extracted_sample_info/snp/{sample}.tsv"
    params:
        position = lambda wildcards: config['important_snpmania_positions'],
        experiment = lambda wildcards: samples['experiment'][wildcards.sample],
        tumour = lambda wildcards: samples['tissue'][wildcards.sample]
    log:
        "logs/reports/extracted_sample_info/snp/{sample}.tsv"
    run:
        extract_snp_file(
            input[0],
            output[0],
            [position_filter(params.position)],
            wildcards.sample,
            params.experiment,
            params.tumour)

def  get_samples():
    return [sample.Index for sample in samples.itertuples()]

rule combine_snp:
    input:
        expand("reports/extracted_sample_info/snp/{sample}.tsv", sample=get_samples())
    output:
        "reports/all.snp.tsv"
    run:
        shell("""cat {input} | awk 'BEGIN{{print("#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tId\tRef_amp\tVar_amp")}} {{if($1!~/#Run/){{print($0)}}}}' > {output}""")
