# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import re

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

# Dictionaries used to extract the correct column in a jSNPmania file
_column_converter_snv = {
    'depth': 0,
    'allele_ratio': 2,
    'nc_number': 3,
    'position': 4,
    'reference': 5,
    'alleles_depth': 6,
    'tumor_A': 10,
    'tumor_G': 11,
    'tumor_C': 12,
    'tumor_T': 13,
    'amplicon_information': -1}

_column_converter_indel = {
    'depth': 0,
    'nc_number': 1,
    'position': 2,
    'allele_ratio': 4,
    'indel': 5,
    'tumor_ind': 6,
    'amplicon_information': -1}

# Order of bases in a jSNPmania file.
_ref_position = {'A': 0, 'G': 1, 'C': 2, 'T': 3}

def convert_jsnpmania_to_annvovar_output(sample_name, output, jsnpmania_variants, jsnpmania_insertion, jsnpmania_deletions, path_nc_to_chr, min_allele_ratio, min_read_depth, amplicon_min_depth):
    chr_to_nc = dict()
    with open(path_nc_to_chr,'r') as nc_to_chr_file:
        for line in nc_to_chr_file:
            if not line.startswith("#"):
                columns = line.split("\t")
                chr_to_nc[columns[0]] = columns[1]
    process_jsnpmania_files(sample_name, output, jsnpmania_variants, jsnpmania_insertion, jsnpmania_deletions, chr_to_nc, min_allele_ratio, min_read_depth, amplicon_min_depth)

def process_jsnpmania_files(sample_name, output, jsnpmania_variants, jsnpmania_insertion, jsnpmania_deletions, nc_to_chr, min_allele_ratio, min_read_depth, amplicon_min_depth):
    """
        Main function, used to process jsnpmania files to generated annovor input files.


        >>> nc_to_chr = {'NC_000001.10': '1', 'NC_000002.11': '2', 'NC_000003.11': '3', 'NC_000004.11': '4', 'NC_000007.13': '7', 'NC_000010.10': '10', 'NC_000012.11': '12', 'NC_000014.8': '14', 'NC_000015.9': '15', 'NC_000017.10': '17', 'NC_000019.9': '19'}
        >>> process_jsnpmania_files("17-1891", "/snakemake-workflows/temp/17-1891.test", "/snakemake-workflows/temp/17-1891.ampliconmapped.variations", "/snakemake-workflows/temp/17-1891.ampliconmapped.insertions", "/snakemake-workflows/temp/17-1891.ampliconmapped.deletions",nc_to_chr,0.01,20,5)
    """
    with open(output,'w') as output_file:
        output_file.write("#Min_tumor_read_depth=" + str(min_read_depth))
        output_file.write("\n#Min_tumor_variant_allele_ratio=" + str(min_allele_ratio))
        output_file.write("\n#Min_amplicon_read_depth=" + str(amplicon_min_depth))
        # Reference dictionary, format: {'NC_000001.10': {'162722758': {'reference': '..', ...}}}
        reference = dict()
        variants = dict()
        # Populate reference dictionary with chromesome identifiers and empty dictionaries
        for keys in nc_to_chr:
            reference[keys] = {}
        # Process jSNPmania SNVs, this will both generate a reference dictionary and a dictionary containing SNV
        with open(jsnpmania_variants,'r') as j_variants:
            # Function used to filter data
            def major_allele_filter(allele_depth, columns): return (
                int(columns[0]) >= min_read_depth and
                int(allele_depth)/int(columns[0]) >= min_allele_ratio and
                (1.0 - float(columns[_column_converter_snv['allele_ratio']])) >= min_allele_ratio)

            for line in j_variants:
                line = line.rstrip()
                if not line.startswith("#") and not line.startswith("\n") != 0:
                    # Populate referebce dictionary
                    reference = extract_ref_variant_info(reference, line, amplicon_min_depth)
                    # Populate variant dictionary
                    alleles = extract_major_snv_allele(sample_name, reference, line, nc_to_chr, amplicon_min_depth, major_allele_filter)
                    for key, info in alleles.items():
                        (nc, start, stop, ref, var) = key.split("#")
                        try:
                            if int(start) in variants[nc]:
                                # Additional variant found for the current position
                                variants[nc][int(start)][key] = info
                            else:
                                # First variant found at the current position
                                variants[nc][int(start)] = {key: info}
                        except KeyError:
                            # First variant found for the current chromosome
                            variants[nc] = {int(start): {key: info}}
        # Print all found SNVs
        for nc in variants:
            for pos in variants[nc]:
                for key, info in variants[nc][pos].items():
                    output_file.write("\n" + info)

        # Iterate through all SNVs and combine the SNVs laying next to each other into Multi Nucleotide Variants (MNVs)
        mb_variants = generate_multi_bp_variants(sample_name,nc_to_chr,variants,amplicon_min_depth)
        # Print all found MNVs
        for key, value in mb_variants.items():
            output_file.write("\n" + value)

        # Process jSNPmania insertions.
        variants = dict()
        with open(jsnpmania_insertion,'r') as j_insertions:
            def insertion_filter(variant_depth, columns): return int(columns[_column_converter_indel['depth']]) >= min_read_depth and int(variant_depth)/int(columns[_column_converter_indel['depth']]) >= min_allele_ratio
            for line in j_insertions:
                line = line.rstrip()
                if not line.startswith("#") and not line.startswith("\n") != 0:
                    alleles = extract_insertion(sample_name, reference, line, nc_to_chr, amplicon_min_depth, insertion_filter)
                    for key, info in alleles.items():
                        (nc, start, stop, ref, var) = key.split("#")
                        try:
                            if int(start) in variants[nc]:
                                # Additional variant found for the current position
                                variants[nc][int(start)][key] = info
                            else:
                                # First variant found at the current position
                                variants[nc][int(start)] = {key: info}
                        except KeyError:
                            # First variant found for the current chromosome
                            variants[nc] = {int(start): {key: info}}
        # Print all found insertions
        for nc in variants:
            for pos in variants[nc]:
                for key, info in variants[nc][pos].items():
                    output_file.write("\n" + info)

        # Process jSNPmania insertions.
        variants = dict()
        with open(jsnpmania_deletions,'r') as j_deletions:
            def deletions_filter(variant_depth, columns): return int(columns[_column_converter_indel['depth']]) >= min_read_depth and int(variant_depth)/int(columns[_column_converter_indel['depth']]) >= min_allele_ratio
            for line in j_deletions:
                line = line.rstrip()
                if not line.startswith("#") and not line.startswith("\n") != 0:
                    alleles = extract_deletions(sample_name, reference, line, nc_to_chr, amplicon_min_depth, deletions_filter)
                    for key, info in alleles.items():
                        (nc, start, stop, ref, var) = key.split("#")
                        try:
                            if int(start) in variants[nc]:
                                if key in variants[nc][int(start)]:
                                    # Additional variant found for the current position
                                    variants[nc][int(start)] = _update_deletion_data(variants[nc][int(start)],{key: info})
                                else:
                                    # First variant found at the current position
                                    variants[nc][int(start)][key] = info
                            else:
                                variants[nc][int(start)] = {key: info}
                        except KeyError:
                            #First variant found for the current chromosome
                            variants[nc] = {int(start): {key: info}}
        # Print all found deletions
        for nc in variants:
            for pos in variants[nc]:
                for key, info in variants[nc][pos].items():
                    output_file.write("\n" + info)


def extract_ref_variant_info(ref, line, amplicon_min_depth=0):
    """
        Extract reference variant information from a jSNPmania variations file.

        :param ref: (dict) that will be populated with reference information, must already
                contain keys for the expected NC numbers, ex {'NC_000001.10': {}, 'NC_000002.8': {}}
        :param line: (string) entry from jSNPmania variations file
        :param amplicon_min_depth: (int) minum read depth for an amplicon if's going to be
                included in result, if set to null no amplicon information will be
                extracted.

        :return
            ref= {
                NC_000001.10: {
                    162722758: {
                        'reference': 'A',
                        'depth': 365,
                       'amp+': 5,
                      'amp-': 1,
                      'ampInfo': 'chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195...'
                    }
                }
            }
    """
    columns = re.split("\t", line.rstrip())
    try:
        nc_number = columns[_column_converter_snv['nc_number']]
        position = columns[_column_converter_snv['position']]
        ref_index = _ref_position[columns[_column_converter_snv['reference']]]
        ref[nc_number][position] = { 'reference': columns[_column_converter_snv['reference']],
                'depth': _extract_allele_counts_snv(columns)[ref_index]}
        if amplicon_min_depth:
            amp_info, amp_p, amp_m = _get_ref_amplicon_information(columns,ref_index, amplicon_min_depth)
            ref[nc_number][position]['amp+'] = amp_p
            ref[nc_number][position]['amp-'] = amp_m
            ref[nc_number][position]['ampInfo'] = amp_info
    except KeyError:
        logger.info("Unexpected reference base detected:" + columns[_column_converter_snv['reference']] + "!")
    return ref


def extract_major_snv_allele(sample, ref, line, nc_to_chr, amplicon_min_depth=0, filter_variants=None):
    """
        Extract major allele information from a JSNPmania variation line, if multiple variants have the same
        allele ratio all of them will be returned.

        :param sample: name of sample
        :param ref: (dict) with reference information, generated by extract_ref_variant_info
        :param line: line from a jSNPmania variants file
        :param nc_to_chr: dict used to convert chr identifier to the required output format
        :param amplicon_min_depth; minimum required number of reads mapped to an amplicon to count it as a valid amplicon, no
                amplicon information will be checked or added if amplicon_min_depth is set to 0 (which is default).
        :param filter_variants: function used to filter variant data

        :return
            {'CHR#START#STOP#REF#VAR': Annovar_input_line,
              ...
            }

    """
    columns = re.split("\t", line.rstrip())

    def var_data(field): return columns[_column_converter_snv[field]]
    alleles = {}
    try:
        ref_data = ref[var_data('nc_number')][var_data('position')]
        variants = _ref_position.copy()
        del variants[var_data('reference')]
        nc_number = var_data('nc_number')
        position = var_data('position')
        for key, allele in variants.items():
            allele_depth = _extract_allele_counts_snv(columns)[allele]
            major_allele_amplicon_information = var_data('amplicon_information').split("|")[allele]
            if filter_variants is None or filter_variants(allele_depth, columns):
                major_alle_info = "\t".join([
                    nc_to_chr[nc_number],
                    position,
                    position,
                    columns[_column_converter_snv['reference']],
                    key,
                    "comments: sample=" + sample +
                    " variantAlleleRatio=" + ("%.16f" % (int(allele_depth)/(int(var_data('depth'))))).rstrip('0') +
                    " alleleFreq=" + str(ref_data['depth']) + "," + str(allele_depth) +
                    " readDepth=" + var_data('depth') +
                    " Tumor_A=" + var_data('tumor_A') +
                    " Tumor_G=" + var_data('tumor_G') +
                    " Tumor_C=" + var_data('tumor_C') +
                    " Tumor_T=" + var_data('tumor_T')])

                if amplicon_min_depth:
                    try:
                        amp_p, amp_m = _count_amplicons(major_allele_amplicon_information, amplicon_min_depth)
                    except ValueError:
                        amp_p, amp_m = 0, 0
                    major_alle_info += " Tumor_var_plusAmplicons=" + str(amp_p) + \
                        " Tumor_var_minusAmplicons=" + str(amp_m) + \
                        " Tumor_ref_plusAmplicons=" + str(ref_data['amp+']) + \
                        " Tumor_ref_minusAmplicons=" + str(ref_data['amp-']) + \
                        " Tumor_var_ampliconInfo=" + major_allele_amplicon_information + \
                        " Tumor_ref_ampliconInfo=" + ref_data['ampInfo']
                alleles[nc_number + '#' + position + "#" + position + "#" + columns[_column_converter_snv['reference']] + "#" + key] = major_alle_info
    except AttributeError:
        logger.info("Unexpected reference base detected" + var_data('reference') + " when extracting major allele!")
        exit(2)
    return alleles


def _get_major_vaf(variants_at_position, sep=r"\s|=", allele_ratio_column=10):
    """
        Extracts the major allele from a dict containing multiple variants.

        :param variants_at_position: a dict {key: data1, key2: data2}, generated using _extract_major_snv_allele
        :param sep: data separator characters(s)
        :param allele_ratio_column: column position where allele ratio can be found

        :return:
            (key: best_allele_ratio_value)
    """
    variant_iterator = iter(variants_at_position.keys())
    key = next(variant_iterator)
    import re
    vaf = float(re.split(sep,variants_at_position[key])[allele_ratio_column] + "\n")
    #ToDo write a test case that goes through this for loop.
    for next_key in variant_iterator:
        next_vaf = float(re.split(r"\s|=",variants_at_position[next_key])[allele_ratio_column])
        if vaf < next_vaf:
            vaf = next_vaf
            key = next_key
    return (key, vaf)


def generate_multi_bp_variants(sample, nc_to_chr, variants, amplicon_min_depth=0):
    """
        Finds variants laying next to each other and passes them to create_multi_bp_variant to generate a multi variant base

        Args:
            sample: name of sample
            nc_to_chr: dict used to convert NC id to chr id
            variants: a dict with all found variants
            amplicon_min_depth:

        Returns:

    """
    multi_bp_variants = dict()
    multi_bp_variant = dict()
    for chr_nc in variants:
        positions = sorted(variants[chr_nc].keys())
        for i in range(0,len(positions)):
            position = positions[i]
            multi_bp_variant[position] = variants[chr_nc][position]
            for j in range(i+1,len(positions)):
                if position + 1 == positions[j]:
                    position = positions[j]
                    multi_bp_variant[position] = variants[chr_nc][position]
                    if len(multi_bp_variant.keys()) > 1:
                        (key, variant) = _create_multi_bp_variant(sample,nc_to_chr[chr_nc],multi_bp_variant,amplicon_min_depth)
                        multi_bp_variants[key] = variant
                else:
                    break

            multi_bp_variant = dict()
    return multi_bp_variants


def _create_multi_bp_variant(sample, chr_name, variants, amplicon_min_depth=0):
    """
        Create Multi Nucleotide Variants (MNV), by combining single snvs laying next to each other, the lowest allele ratio found
        in a snv building up the

        Args:
            sample: name of sample
            chr_name: chromosome id/name
            variants: collection of variants laying next to each other
            amplicon_min_depth: (int) min required read depth for an amplicon to be included

        Returns:
            ('CHR#START#STOP#REF#VAR', Annovar_input_line,)
    """
    positions = iter(sorted(variants.keys()))
    start_position = next(positions)
    (key, best_vaf) = _get_major_vaf(variants[start_position])
    (chr_nc, start, stop, new_ref, new_var ) = key.split("#")
    info = variants[start_position][key]
    first_info = info

    for last_position in positions:
        #Sorting is done to be able to emulate sera result
        from collections import OrderedDict
        (key, vaf) = _get_major_vaf(OrderedDict(sorted(variants[last_position].items(), key=lambda t: t[0], reverse=True)))
        (chr_nc, start, stop, next_ref, next_var ) = key.split("#")
        new_ref +=  next_ref
        new_var += next_var
        if vaf < best_vaf:
            best_vaf = vaf
            info = variants[last_position][key]
    import re
    columns = re.split(r"\s|=",info)
    columns_incorrect = re.split(r"\s|=",first_info)
    mb_variant = "\t".join([chr_name,str(start_position),str(last_position),new_ref,new_var,
            "comments: sample=" + sample + " variantAlleleRatio=" +
              str(best_vaf) + " alleleFreq=" + columns[11] + " readDepth=" + columns[13] +
              " Tumor_A=- Tumor_G=- Tumor_C=- Tumor_T=-"])
    if amplicon_min_depth:
        # ToDo change to columns to fix incorrect printing of amplicon info,
        # the current SERA has a bug
        mb_variant += " Tumor_var_plusAmplicons=" + columns_incorrect[23] + \
                " Tumor_var_minusAmplicons=" + columns_incorrect[25] + \
                " Tumor_ref_plusAmplicons=" + columns_incorrect[27] + \
                " Tumor_ref_minusAmplicons=" + columns_incorrect[29] + \
                " Tumor_var_ampliconInfo=" + columns_incorrect[31] + \
                " Tumor_ref_ampliconInfo=" + columns_incorrect[33]
    return (chr_nc + "#" + str(start_position) + "#" + str(last_position) + "#" + new_ref + "#" + new_var, mb_variant)


def _extract_allele_counts_snv(columns):
    return [int(allele_d) for allele_d in columns[_column_converter_snv['alleles_depth']].split("|")]


def extract_deletions(sample, ref, line, nc_to_chr, amplicon_min_depth=0, deletion_filter=None):
    """
        Extract deletions information from a JSNPmania deletion line.

        :param sample: (String) sample name
        :param ref: (dict) with reference information, generated by extract_ref_variant_info
        :param line: line from a jSNPmania deletions file
        :param nc_to_chr: dict used to convert chr identifier to the required output format
        :param amplicon_min_depth; minimum required number of reads mapped to an amplicon to count it as a valid amplicon, no
                amplicon information will be checked or added if amplicon_min_depth is set to 0 (which is default).
        :param deletion_filter: function used to filter deletion data

        :return
            {'CHR#START#STOP##-': Annovar_input_line,
              ...
            }
    """
    columns = re.split(r"\t", line.rstrip())

    def var_data(field): return columns[_column_converter_indel[field]]
    ref_data = ref[var_data('nc_number')][var_data('position')]
    position = int(var_data('position'))
    nc_number = var_data('nc_number')
    deletions = re.split(r"\|", var_data('indel'))
    extracted_deletions = {}
    for deletion in deletions:
        if deletion != '0':
            depth, variant = re.search(r"^(\d+)\(([0-9,-]+)\)$", deletion).groups()
            if deletion_filter is None or deletion_filter(depth ,columns):
                from_bp, to_bp = map(int,variant.split(","))
                deleted_bases = "".join([ref[nc_number][str(position + p)]['reference'] for p in range(from_bp,to_bp + 1)])
                allele_ratio = int(depth)/int(var_data('depth'))
                key = nc_number + "#" + str(from_bp + position) + "#" + str(to_bp + position) + "##-"
                deletion_info = "\t".join([
                    nc_to_chr[nc_number],
                    str(from_bp + position), str(to_bp + position),
                    deleted_bases,
                    '-',
                    "comments: sample=" + sample +
                    " variantAlleleRatio=" + str(allele_ratio) +
                    " alleleFreq=" + str(ref_data['depth']) + "," + depth +
                    " readDepth=" + var_data('depth') +
                    " Tumor_Del=" + var_data('tumor_ind')])
                ampliconInfo = None
                if amplicon_min_depth:
                    for amp_info in var_data('amplicon_information').split("|"):
                        if amp_info.startswith("(" + variant + ")"):
                            amp_p, amp_m = _count_amplicons(amp_info,amplicon_min_depth)
                            ampliconInfo = " Tumor_var_plusAmplicons=" + str(amp_p) + \
                                           " Tumor_var_minusAmplicons=" + str(amp_m) + \
                                           " Tumor_ref_plusAmplicons=" + str(ref_data['amp+']) + \
                                           " Tumor_ref_minusAmplicons=" + str(ref_data['amp-']) + \
                                           " Tumor_var_ampliconInfo=" + amp_info + \
                                           " Tumor_ref_ampliconInfo=" + ref_data['ampInfo']
                            break;
                    if ampliconInfo:
                        deletion_info += ampliconInfo
                    else:
                        deletion_info += " Tumor_var_plusAmplicons=-" + \
                                         " Tumor_var_minusAmplicons=-" + \
                                         " Tumor_ref_plusAmplicons=" + str(ref_data['amp+']) + \
                                         " Tumor_ref_minusAmplicons=" + str(ref_data['amp-']) + \
                                         " Tumor_var_ampliconInfo=-" + \
                                         " Tumor_ref_ampliconInfo=" + ref_data['ampInfo']
                extracted_deletions[key] = deletion_info
    return extracted_deletions


def _update_deletion_data(variants, deletions):
    """
        Compare new deletion with previously found deletions and update variants with the new deletion if
        has a better allele ratio

        :param variants: (dict) previously found deletions
        :param deletions: new deletion

        :return:
    """
    def extract_info(info):
        return dict(tuple(comment.split("=")) for comment in info.split("\t")[-1].split(" ")
                    if not comment.startswith("comment"))
    for key, info in deletions.items():
        if key in variants:
            existing_variant_info = extract_info(variants[key])
            new_variant_info = extract_info(info)
            if float(new_variant_info['variantAlleleRatio']) > float(existing_variant_info['variantAlleleRatio']):
                variants[key] = info
        else:
            variants[key] = info
    return variants


def extract_insertion(sample, ref, line, nc_to_chr, amplicon_min_depth=0, insertion_filter=None):
    """
        Extract insertion information from a JSNPmania insertion line.

        :param sample: (String) sample name
        :param ref: (dict) with reference information, generated by extract_ref_variant_info
        :param line: line from a jSNPmania insertion file
        :param nc_to_chr: dict used to convert chr identifier to the required output format
        :param amplicon_min_depth; minimum required number of reads mapped to an amplicon to count it as a valid
                amplicon, no amplicon information will be checked or added if amplicon_min_depth is set to 0 (which is
                default).
        :param deletion_filter: function used to filter deletion data

        :return
            {'CHR#START#STOP##-': Annovar_input_line,
              ...
            }
    """
    columns = re.split(r"\t", line.rstrip())

    def var_data(field): return columns[_column_converter_indel[field]]
    ref_data = ref[var_data('nc_number')][var_data('position')]
    insertions = re.split(r"\|", var_data('indel'))
    position = var_data('position')
    nc_number = var_data('nc_number')
    extracted_insertions = {}
    ref_depth = ref[nc_number][position]['depth']
    for index, insertion in enumerate(insertions):
        if insertion != '0':
            depth, variant = re.search(r"^(\d+)([A-Z]+$)", insertion).groups()
            if insertion_filter is None or insertion_filter(depth,columns):
                allele_ratio = int(depth)/int(var_data('depth'))
                insertion_info = "\t".join([
                    nc_to_chr[nc_number],
                    position, position,
                    '-',
                    variant,
                    "comments: sample=" + sample +
                    " variantAlleleRatio=" + str(allele_ratio) +
                    " alleleFreq=" + str(ref_depth) + "," + depth +
                    " readDepth=" + var_data('depth') +
                    " Tumor_Ins=" + var_data('tumor_ind')])
                if amplicon_min_depth:
                    for amp_info in var_data('amplicon_information').split("|"):
                        if amp_info.startswith(variant + ":"):
                            amp_p, amp_m = _count_amplicons(amp_info,amplicon_min_depth)
                            insertion_info += " Tumor_var_plusAmplicons=" + str(amp_p) + \
                                              " Tumor_var_minusAmplicons=" + str(amp_m) + \
                                              " Tumor_ref_plusAmplicons=" + str(ref_data['amp+']) + \
                                              " Tumor_ref_minusAmplicons=" + str(ref_data['amp-']) + \
                                              " Tumor_var_ampliconInfo=" + amp_info + \
                                              " Tumor_ref_ampliconInfo=" + ref_data['ampInfo']
                            break
                extracted_insertions[nc_number + "#" + position + "#" + position + "#-#" + variant] = insertion_info
    return extracted_insertions


def _count_amplicons(amplicon_information, min_depth=0):
    """
        Count number of HaloPlex amplicons for each strand that have fulfilled
        the minimum read depth.

        :param amplicon_information: (string) a jSNPmania generate amplicon information string.
                SNV example: chr1:162722740-162722971:+:1|chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195
                Insertion/Deletion example: (0,0):chr1:162722740-162722971:+:1
        :param min_depth: (int) min required read depth for an amplicon to be included

        :return
            (int,int) --> (number of plus strand amplicons, number of minus strand amplicons)
    """
    plus_strand_amplicons, minus_strand_amplicons = 0, 0
    for amplicon in amplicon_information.split("#"):
        strand, depth = amplicon.split(":")[-2:]
        if int(depth) >= min_depth:
            if strand == "+":
                plus_strand_amplicons+=1
            else:
                minus_strand_amplicons+=1
    return (plus_strand_amplicons, minus_strand_amplicons)


def _get_ref_amplicon_information(columns, ref_index, amplicon_min_depth):
    """
        Extract amplicon information from jSNPmania variation file line.

        :param columns: list of data values, originating from a jSNPmania variations file line.
        :param ref_index: ref index positon, i.e 0,1,2 or 3
        :param amplicon_min_depth:: min allowed depth for

        :return: ("amplicon information", number of amplicon for plus strand, number of amplicon for minus strand)
    """
    try:
        amplicon_information = columns[_column_converter_snv['amplicon_information']].split("|")[ref_index]
        amp_p, amp_m = _count_amplicons(amplicon_information, amplicon_min_depth)
        return amplicon_information, amp_p, amp_m
    except ValueError:
        return '0', 0, 0
