# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


_umi_r2_head_input_bam = "mapped/{sample}.sorted.bam"
try:
    _umi_r2_head_input_bam = umi_r2_head_input_bam
except:
    pass

_umi_r2_head_input_umi = "mapped/{sample}.sorted.bam"
try:
  _umi_r2_head_input_umi = umi_r2_head_input_umi
except:
    pass

_umi_r2_head_output = "mapped/{sample}.amplicon_annotated.sorted.bam
try:
    _umi_r2_head_output = umi_r2_head_output
except:
    pass

rule annotate_bam_with_umis_part:
    input:
        bam=_umi_r2_head_input_bam,
        umi=_umi_r2_head_input_umi
    output:
        bam=_umi_r2_head_output
    log:
        "logs/umis/{sample}.{unit}.{part}.fgbioAnnoBam.txt"
    wrapper:
        "fgbio/bio/fgbio/annotatebamwithumis"
