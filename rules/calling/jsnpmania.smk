# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.lib.common.utils import get_bam_file

_jsnpmania_input = lambda wildcards: "mapped/" + get_bam_file(wildcards, samples)
try:
    _jsnpmania_input = jsnpmania_input
except:
    pass

_jsnpmania_variations_output = "jsnpmania/{sample}.variations"
try:
    _jsnpmania_variations_output = jsnpmania_variations_output
except:
    pass

_jsnpmania_insertions_output = "jsnpmania/{sample}.insertions"
try:
    _jsnpmania_insertions_output = jsnpmania_insertions_output
except:
    pass

_jsnpmania_deletions_output = "jsnpmania/{sample}.deletions"
try:
    _jsnpmania_deletions_output = jsnpmania_deletions_output
except:
    pass

rule run_jsnpmania:
    input:
        _jsnpmania_input
    output:
        variations = _jsnpmania_variations_output,
        insertions = _jsnpmania_insertions_output,
        deletions = _jsnpmania_deletions_output
    params:
        path_jsnpmania = config.get('path_jsnpmania',"jsnpmania.jar"),
        path_jsnpmania_header = config['path_jsnpmania_header'],
        ref_file = lambda wildcards: samples['snpseq_file'][wildcards.sample],
        flags = config["jsnpmania_flags"]
    log:
        "logs/jsnpmania/{sample}.jsnpmania.log"
    wrapper:
        "master/bio/jsnpmania"#"master/bio/jsnpmania"
