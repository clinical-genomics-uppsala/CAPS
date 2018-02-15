# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


def convert_jsnpseq_snv_entry(ref, line,amplicon_mapped=True):
    column_converter = {'read_depth': 0, 'ref_read_depth': 1, 'allele_ratio': 2
                        'nc_number': 3, 'position': 4, 'reference': 5,
                        'alleles_depth': 6, 'amplicon_information': 14}
    ref_position = {'A': 0, 'G': 1, 'C': 2, 'T': 3)
    columns = line.split("\t")
    alleles = columns.split("\")
    try:
        nc_number = columns[column_converter['nc_number']]]
        position = columns[column_converter['position']]]
        ref_index = ref_position[columns[column_converter['reference']]]]
        ref[nc_number][position] = columns[column_converter['reference']]]
        ref[nc_number][position] = columns[column_converter['alleles_depth']]].split("|")[ref_index]
        if amplicon_mapped:

def count_amplicons(amplicon_information, type, min_amplicon_depth=0):
    plus_strand_amplicons, minus_strand_amplicons = 0, 0
    for amplicon in amplicon_information.split("#"):
        strand, depth = amplicon.split(":")[-2:]
        if depth >= min_amplicon_depth:
            if strand == "+":
                plus_strand_amplicons++
            else:
                minus_strand_amplicons++
