from scripts.lib.common.data.report.wp1 import calc_amplication

rule calc_amplication:
    input:
        "jsnpmania/{sample}.variations"
    output:
        "reports/amplification/{sample}.amplification.tsv"
    params:
        background_regions = lambda wildcards: samples['background_file'][wildcards.sample],
        amplication_regions = lambda wildcards: samples['amplification_file'][wildcards.sample],
        tissue = lambda wildcards: samples['tissue'][wildcards.sample],
        chr_to_nc = config["chr_to_nc_file"]
    run:
        calc_amplication(
                  wildcards.sample,
                  params.tissue,
                  output[0],
                  input[0],
                  params.amplication_regions,
                  params.background_regions,
                  params.chr_to_nc
                )
