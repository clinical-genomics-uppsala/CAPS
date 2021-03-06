# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

_annovar_input = "annovar/{sample}.{type}.annovarInput"
try:
    _annovar_input = annovar_input
except:
    pass

_annovar_output = "annovar_output/{sample}.{type}.singleSample.annovarOutput"
try:
    _annovar_output = annovar_output
except:
    pass

rule run_annovar_table:
    input:
         _annovar_input
    output:
        temp("annovar/{sample}.{type}.annovarInput.hg19_multianno.txt")
    params:
        table_annovar = config['annovar']["bin"],
        database = config['annovar']["database"],
        protocols = ",".join([protocol for protocol in config['annovar']['protocols']]),
        operations = ",".join([config['annovar']['protocols'][protocol]['operation'] for protocol in config['annovar']['protocols']]),
        arguments = ",".join([config['annovar']['protocols'][protocol]['argument'] for protocol in config['annovar']['protocols']]),
        build_version = config['annovar']['build_version']
    log:
        "logs/annovar/{sample}.{type}.table_annovar.log"
    run:
        shell("perl {params.table_annovar} {input} {params.database} -protocol {params.protocols} -operation {params.operations} -nastring \"-\" -otherinfo -buildver {params.build_version} -remove -arg {params.arguments} 2> {log}")

from scripts.lib.common.data.parser.annovar import process_annovar_multianno_file

def field_converter(value):
    if value == "yes":
        return True
    elif value == "no":
         False
    else:
        raise Exception("Unhandled field value: " + value)

def get_sample_value(field, sample, converter, default):
    if field in samples and sample in samples[field]:
        return converter(samples[field][sample])
    else:
        return default

rule create_annovar_output:
    input:
         "annovar/{sample}.{type}.annovarInput.hg19_multianno.txt"
    output:
        _annovar_output
    params:
        tumor_normal = lambda wildcards: get_sample_value('tumor_normal', wildcards.sample, converter=field_converter),
        amplicon_mapped = lambda wildcards: get_sample_value('amplicon_mapped', wildcards.sample, converter=field_converter)
    log:
        "logs/annovar/{sample}.{type}.table_annovar.log"
    run:
        process_annovar_multianno_file(
            output[0],
            input[0],
            params.tumor_normal,
            params.amplicon_mapped)

#rule run_annovar_pindel_table:
#    input:
#         "pindel_annovar/{sample}.pindel.filtered.annovarInput"
#    output:
#        "pindel_annovar/{sample}.pindel.filtered.annovarInput.hg19_multianno.txt"
#    params:
#        table_annovar = config['annovar']["bin"],
#        database = config['annovar']["database"],
#        protocols = ",".join([protocol for protocol in config['annovar']['protocols']]),
#        operations = ",".join([config['annovar']['protocols'][protocol]['operation'] for protocol in config['annovar']['protocols']]),
#        arguments = ",".join([config['annovar']['protocols'][protocol]['argument'] for protocol in config['annovar']['protocols']]),
#        build_version = config['annovar']['build_version']
#    log:
#        "logs/pindel_annovar/{sample}.table_annovar.log"
#    run:
#      shell("perl {params.table_annovar} {input} {params.database} -protocol {params.protocols} -operation {params.operations} -nastring \"-\" -otherinfo -buildver {params.build_version} -remove -arg {params.arguments} 2> {log}")
#
#rule create_pindel_annovar_output:
#    input:
#         "pindel_annovar/{sample}.pindel.filtered.annovarInput.hg19_multianno.txt"
#    output:
#        "pindel_annovar_output/{sample}.pindel.singleSample.annovarOutput"
#    params:
#        tumor_normal = lambda wildcards: get_sample_value('tumor_normal', wildcards.sample, converter=field_converter),
#        amplicon_mapped = lambda wildcards: get_sample_value('amplicon_mapped', wildcards.sample, converter=field_converter)
#    log:
#        "logs/pindel_annovar/{sample}.table_annovar.log"
#    run:
#        process_annovar_multianno_file(
#            output[0],
#            input[0],
#            params.tumor_normal,
#            params.amplicon_mapped)
