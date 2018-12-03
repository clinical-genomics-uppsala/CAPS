# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def collapse_reads_input():
    try:
      return collapse_reads_input
    except:
      return "mapped/{sample}.{unit}.{part}.primerclip.umi.bam"

def collapse_reads_output():
    try:
      return collapse_reads_output
    except:
      return "mapped/consensus/{sample}.consensus.bam"

rule umi_cr_revertsam:
    input:
        collapse_reads_input()
    output:
        temp("mapped/consensus/{sample}.{unit}.{part}.sanitised.bam")
    params:
        extra="SANITIZE=true REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false"
    wrapper:
        "picard-revertsam/bio/picard/revertsam"

rule set_mq_information_fgbio:
    input:
        bam="mapped/consensus/{sample}.{unit}.{part}.sanitised.bam"
    output:
        bam=temp("mapped/consensus/{sample}.{unit}.{part}.sanitised.setmateinfo.bam")
    log:
        "logs/umis/collapse/{sample}.{unit}.{part}.fgbioSetMQInfo.txt"
    wrapper:
        "fgbio/bio/fgbio/setmateinformation"

def get_bam_files(units, config):
  num_splits = config.get("cgu_accel_num_fastq_split", config.get("num_fastq_split", 1))
  if num_splits > 1:
    return [ unit + ".%04d" % part for part in range(0,num_splits) for unit in units]
  else:
    return units

rule merge_bam_files_fgbio:
    input:
        lambda wildcards: expand("mapped/consensus/{sample_unit}.sanitised.setmateinfo.bam", sample_unit=get_bam_files(get_units(wildcards,units),config))
    output:
        temp("mapped/consensus/{sample}.fgbio.merged.bam")
    threads: 8
    wrapper:
        "fgbio/bio/samtools/merge"

rule sort_sanitised_queryname_fgbio:
    input:
      "mapped/consensus/{sample}.fgbio.merged.bam"
    output:
        temp("mapped/consensus/{sample}.fgbio.merged.setmateinfo.qsorted.bam")
    log:
        "logs/umi/collapse/qsort/{sample}.log"
    params:
        sort_order="queryname",
        extra="VALIDATION_STRINGENCY=LENIENT -Xms500m -Xmx20g"
    threads: 8
    wrapper:
        "master/bio/picard/sortsam"

rule set_groupreads_by_umi_fgbio:
    input:
        "mapped/consensus/{sample}.fgbio.merged.setmateinfo.qsorted.bam"
    output:
        temp("mapped/consensus/{sample}.merged.sanitised.setmateinfo.groupreads.bam")
    log:
        "logs/umis/collapse/{sample}.merged.fgbioGroup.txt"
    params:
        extra="-s adjacency --edits 1"
    wrapper:
        "fgbio/bio/fgbio/groupreadsbyumi"

rule create_consensus_reads_m3:
    input:
        "mapped/consensus/{sample}.merged.sanitised.setmateinfo.groupreads.bam"
    output:
        temp("mapped/consensus/{sample}.merged.sanitised.setmateinfo.groupread.m3s.bam")
    log:
       "logs/umis/collapse/{sample}.merged.fgbioCMCR-M3.txt"
    params:
        extra="-M 3"
    wrapper:
        "fgbio/bio/fgbio/callmolecularconsensusreads"

rule create_consensus_reads_m1:
    input:
        "mapped/consensus/{sample}.merged.sanitised.setmateinfo.groupreads.bam"
    output:
        temp("mapped/consensus/{sample}.merged.sanitised.setmateinfo.groupread.m1s.bam")
    log:
       "logs/umis/collapse/{sample}.merged.fgbioCMCR-M1.txt"
    params:
        extra="-M 1"
    wrapper:
        "fgbio/bio/fgbio/callmolecularconsensusreads"

rule create_reads_from_bam_m3:
    input:
        "mapped/consensus/{sample}.merged.sanitised.setmateinfo.groupread.m3s.bam"
    output:
        fastq1=temp("mapped/consensus/fastq/{sample}.consensus.m3.R1.fq"),
        fastq2=temp("mapped/consensus/fastq/{sample}.consensus.m3.R2.fq")
    wrapper:
        "picard-samtofastq/bio/picard/samtofastq"

rule bwa_concensus_m3:
    input:
        reads=["mapped/consensus/fastq/{sample}.consensus.m3.R1.fq", "mapped/consensus/fastq/{sample}.consensus.m3.R2.fq"]
    output:
        temp("mapped/consensus/{sample}.consensus.m3.sam")
    log:
        "logs/bwa_mem/{sample}.consensus.log"
    threads: 3
    params:
        index=config['reference_genome'],
        extra=lambda wildcards: r"-M -R '@RG\tID:" + get_now() + "_" + wildcards.sample + r"\tSM:" + wildcards.sample + r"\tPL:illumina'",
        sort="samtools",
        sort_order="queryname",
        sort_extra="-@ 3"
    wrapper:
        "master/bio/bwa/mem"

rule primerclip_consensus:
    input:
        sam="mapped/consensus/{sample}.consensus.m3.sam",
        master_file="master.bed"
    output:
        sam=temp("mapped/consensus/{sample}.consensus.m3.primerclip.sam")
    log:
        "logs/mapped/primerclip/{sample}.consensus.m3.primerclip.log"
    wrapper:
        "primerclip/bio/primerclip"

rule replace_rg_concensus_collapse:
    input:
        "mapped/consensus/{sample}.consensus.m3.primerclip.sam"
    output:
        collapse_reads_output()
    log:
        "logs/picard/replace_rg/{sample}.log"
    params:
        "RGPL=illumina RGPU={sample} RGSM={sample} RGID=Accel  RGLB=MID SO=coordinate VALIDATION_STRINGENCY=STRICT CREATE_INDEX=TRUE"
    wrapper:
        "master/bio/picard/addorreplacereadgroups"
