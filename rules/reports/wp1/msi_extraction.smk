# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.data.report.wp1 import extract_msi_markers

rule extract_msi_markers:
    input:
        "reports/{sample}.filteredMutations.tsv"
    output:
        "reports/extracted_sample_info/msi/{sample}.msiMarkers.tsv"
    params:
        tgfbr2_1bp = 0.01,
        tgfbr2_2bp = 0.01,
        acvr2a_1bp = 0.01,
        acvr2a_2bp = 0.01
    log:
        "logs/reports/extracted_sample_info/msi//{sample}.msiMarkers.log"
    run:
        extract_msi_markers(
            input[0],
            output[0],
            params.tgfbr2_1bp,
            params.tgfbr2_2bp,
            params.acvr2a_1bp,
            params.acvr2a_2bp)

def extract_colon_samples():
    sample_list = []
    for row in samples.itertuples():
        target_type = samples.get("tissue",{}).get(row.Index)
        if target_type == "colon":
            sample_list.append(row.Index)
    return sample_list

rule combine_msi_markers:
    input:
        expand("reports/extracted_sample_info/msi/{sample}.msiMarkers.tsv", sample=extract_colon_samples())
    output:
        "reports/all.msiMarkers.tsv"
    run:
        if len(input) > 1:
            shell("""cat {input} | awk 'BEGIN{{NUM=0}}{{if(!/^#/){{print($0);}}else if(NUM == 0){{print($0);}}NUM=NUM+1;}}' > {output}""")
        elif len(input) == 1:
            shell("ln -sr {input} {output}")
        else:
            shell("touch {output}")
