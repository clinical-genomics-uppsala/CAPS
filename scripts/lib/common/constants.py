trimming = {
    'cutadapt': {
        'fq1': "R1.cutadapt.fastq.gz",
        'fq2': "R2.cutadapt.fastq.gz"
    },
    'swift': {
        'fq1': "R1.trimmomatic_cutadapt.fastq.gz",
        'fq2': "R2.trimmomatic_cutadapt.fastq.gz"
    },
}

bam_file_ending= {
    'haloplex': "amplicon_annotated.sorted.bam",
    'swift': "sorted.bam",
    'default': "sorted.bam"
}
