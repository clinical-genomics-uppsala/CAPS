# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.data.parser.jsnpmania import extract_egfr, g719_filter, t790m_filter

rule extract_egfr_t790m:
    input:
        "jsnpmania/{sample}.variations"
    output:
        "reports/extracted_sample_info/egfr/{sample}.egfr_t790m.tsv"
    params:
        min_depth = 5,
        tissue = lambda wildcards: samples['tissue'][wildcards.sample],
        experiment = lambda wildcards: samples['experiment'][wildcards.sample]
    log:
        "logs/reports/extracted_sample_info/egfr/{sample}.egfr_t790m.log"
    run:
        extract_egfr(
            input[0],
            output[0],
            params.min_depth,
            wildcards.sample,
            params.tissue,
            params.experiment,
            t790m_filter)

rule extract_egfr_g719:
    input:
        "jsnpmania/{sample}.variations"
    output:
        "reports/extracted_sample_info/egfr/{sample}.egfr_g719.tsv"
    params:
        min_depth = 5,
        tissue = lambda wildcards: samples['tissue'][wildcards.sample],
        experiment = lambda wildcards: samples['experiment'][wildcards.sample]
    log:
        "logs/reports/extracted_sample_info/egfr/{sample}.egfr_g719.log"
    run:
        extract_egfr(
            input[0],
            output[0],
            params.min_depth,
            wildcards.sample,
            params.tissue,
            params.experiment,
            g719_filter)

def extract_lung_samples():
    sample_list = []
    for row in samples.itertuples():
        target_type = samples.get("tissue",{}).get(row.Index)
        if target_type == "lung":
            sample_list.append(row.Index)
    return sample_list

rule combine_egfr_t790m:
    input:
        expand("reports/extracted_sample_info/egfr/{sample}.egfr_t790m.tsv", sample=extract_lung_samples())
    output:
        "reports/EGFR_T790M.allSamples.tsv"
    run:
        if len(input) > 0:
            shell("""cat {input} | awk 'BEGIN{print "#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp"} {if($1!~/#Run/){print $0}}' > {output}""")
        else:
            shell("touch {output}")

rule combine_egfr_g719:
    input:
        expand("reports/extracted_sample_info/egfr/{sample}.egfr_g719.tsv", sample=extract_lung_samples())
    output:
        "reports/EGFR_G719.allSamples.tsv"
    run:
        if len(input) > 0:
            shell("""cat {input} | awk 'BEGIN{print "#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp"} {if($1!~/#Run/){print $0}}' > {output}""")
        else:
            shell("touch {output}")
