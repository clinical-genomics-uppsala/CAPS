# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import re

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

_column_converter_snv = {'depth': 0, 'allele_ratio': 2, 'nc_number': 3, 'position': 4, 'reference': 5,
                    'alleles_depth': 6, 'amplicon_information': -1, 'tumor_A': 10,
                     'tumor_G': 11, 'tumor_C': 12, 'tumor_T': 13}

_column_converter_indel = {'depth': 0, 'nc_number': 1, 'position': 2, 'allele_ratio': 4, 'indel': 5, 'tumor_ind': 6, 'amplicon_information': -1}

_ref_position = {'A': 0, 'G': 1, 'C': 2, 'T': 3}

def process_jsnpmania_files(sample_name, output, jsnpmania_variants, jsnpmania_insertion, jsnpmania_deletions, nc_to_chr, min_allele_ratio, min_read_depth, amplicon_min_depth):
    """
        >>> nc_to_chr = {'NC_000001.10': '1', 'NC_000002.11': '2', 'NC_000003.11': '3', 'NC_000004.11': '4', 'NC_000007.13': '7', 'NC_000010.10': '10', 'NC_000012.11': '12', 'NC_000014.8': '14', 'NC_000015.9': '15', 'NC_000017.10': '17', 'NC_000019.9': '19'}
        >>> process_jsnpmania_files("17-1891", "/snakemake-workflows/temp/17-1891.test", "/snakemake-workflows/temp/17-1891.ampliconmapped.variations", "/snakemake-workflows/temp/17-1891.ampliconmapped.insertions", "/snakemake-workflows/temp/17-1891.ampliconmapped.deletions",nc_to_chr,0.01,20,5)
    """
    with open(output,'w') as output_file:
        output_file.write("#Min_tumor_read_depth=" + str(min_read_depth))
        output_file.write("\n#Min_tumor_variant_allele_ratio=" + str(min_allele_ratio))
        output_file.write("\n#Min_amplicon_read_depth=" + str(amplicon_min_depth))
        reference = dict()
        variants = dict()
        for keys in nc_to_chr:
            reference[keys] = {}
        with open(jsnpmania_variants,'r') as j_variants:
            major_allele_filter = lambda allele_depth, columns:  \
                int(columns[0]) >= 20 and  \
                int(allele_depth)/int(columns[0]) >= 0.01 and \
                (1.0 - float(columns[_column_converter_snv['allele_ratio']])) >= 0.01
            for line in j_variants:
                line = line.rstrip()
                if not line.startswith("#") and not line.startswith("\n") != 0:
                    reference = _extract_ref_variant_info(reference, line, amplicon_min_depth)
                    alleles = _extract_major_snv_allele(sample_name, reference, line, nc_to_chr, amplicon_min_depth, major_allele_filter)
                    for key, info in alleles.items():
                        (nc, pos, pos, ref, var) = key.split("#")
                        #try:
                        #    variants[nc][pos][key] = info
                        try:
                            if int(pos) in variants[nc]:
                                variants[nc][int(pos)][key] = info
                            else:
                                variants[nc][int(pos)] = {key: info}
                        except KeyError:
                            variants[nc] = {int(pos): {key: info}}
        for nc in variants:
            for pos in variants[nc]:
                for key, info in variants[nc][pos].items():
                    output_file.write("\n" + info)
        mb_variants = check_for_multi_bp_variants(sample_name,nc_to_chr,variants,amplicon_min_depth)
        for key, value in mb_variants.items():
            output_file.write("\n" + value)
        variants = dict()
        with open(jsnpmania_insertion,'r') as j_insertions:
            insertion_filter = lambda variant_depth, columns: int(columns[_column_converter_indel['depth']]) >= min_read_depth and int(variant_depth)/int(columns[_column_converter_indel['depth']]) >= min_allele_ratio
            for line in j_insertions:
                line = line.rstrip()
                if not line.startswith("#") and not line.startswith("\n") != 0:
                    alleles = _extract_insertion(sample_name, reference, line, nc_to_chr, amplicon_min_depth, insertion_filter)
                    for key, info in alleles.items():
                        (nc, pos, pos, ref, var) = key.split("#")
                        #try:
                        #    variants[nc][pos][key] = info
                        try:
                            if int(pos) in variants[nc]:
                                variants[nc][int(pos)][key] = info
                            else:
                                variants[nc][int(pos)] = {key: info}
                        except KeyError:
                            variants[nc] = {int(pos): {key: info}}
        for nc in variants:
            for pos in variants[nc]:
                for key, info in variants[nc][pos].items():
                    output_file.write("\n" + info)
        variants = dict()
        with open(jsnpmania_deletions,'r') as j_deletions:
            deletions_filter = lambda variant_depth, columns: int(columns[_column_converter_indel['depth']]) >= min_read_depth and int(variant_depth)/int(columns[_column_converter_indel['depth']]) >= min_allele_ratio
            for line in j_deletions:
                line = line.rstrip()
                if not line.startswith("#") and not line.startswith("\n") != 0:
                    alleles = _extract_deletions(sample_name, reference, line, nc_to_chr, amplicon_min_depth, insertion_filter)
                    for key, info in alleles.items():
                        (nc, pos, pos, ref, var) = key.split("#")
                        #try:
                        #    variants[nc][pos][key] = info
                        try:
                            if int(pos) in variants[nc]:
                                variants[nc][int(pos)][key] = info
                            else:
                                variants[nc][int(pos)] = {key: info}
                        except KeyError:
                            variants[nc] = {int(pos): {key: info}}
        for nc in variants:
            for pos in variants[nc]:
                for key, info in variants[nc][pos].items():
                    output_file.write("\n" + info)

def get_major_vaf(variants_at_position):
    variant_iterator = iter(variants_at_position.keys())
    key = next(variant_iterator)
    import re
    vaf = float(re.split("\s|=",variants_at_position[key])[9])
    for next_key in variant_iterator:
        next_vaf = float(re.split("\s|=",variants_at_position[next_key])[9])
        if vaf < next_vaf:
            vaf = next_vaf
            key = next_key
    return (key, vaf)

def create_multi_bp_variant(sample, chr, variants, amplicon_min_depth=0):
    positions = iter(sorted(variants.keys()))
    start_position = next(positions)
    (key, best_vaf) = get_major_vaf(variants[start_position])

    (chr_nc, start, stop, new_ref, new_var ) = key.split("#")
    info = variants[start_position][key]
    first_info = info
    for last_position in positions:
        (key, vaf) = get_major_vaf(variants[last_position])
        (chr_nc, start, stop, next_ref, next_var ) = key.split("#")
        new_ref +=  next_ref
        new_var += next_var
        if vaf < best_vaf:
            best_vaf = vaf
            info = variants[last_position][key]
    import re
    columns = re.split("\s|=",info)
    columns_incorrect = re.split("\s|=",first_info)
    mb_variant = "\t".join([chr,str(start_position),str(last_position),new_ref,new_var,
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

def check_for_multi_bp_variants(sample, nc_to_chr, variants, amplicon_min_depth=0):
    multi_bp_variants = dict()
    multi_bp_variant = dict()
    for chr_nc in variants:
        positions = sorted(variants[chr_nc].keys())
        for i in range(0,len(positions)):
            position = positions[i]
            multi_bp_variant[position] = variants[chr_nc][position]
            for j in range(i+1,len(positions)):
                if(position + 1 == positions[j]):
                    position = positions[j]
                    multi_bp_variant[position] = variants[chr_nc][position]
                    if len(multi_bp_variant.keys()) > 1:
                        (key, variant) = create_multi_bp_variant(sample,nc_to_chr[chr_nc],multi_bp_variant,amplicon_min_depth)
                        multi_bp_variants[key] = variant
                else:
                    break

            #if len(multi_bp_variant.keys()) > 1:
            #    (key, variant) = create_multi_bp_variant(sample,nc_to_chr[chr_nc],multi_bp_variant,amplicon_min_depth)
            #    multi_bp_variants[key] = variant
            multi_bp_variant = dict()
    return multi_bp_variants


def _extract_major_snv_allele(sample, ref, line, nc_to_chr, amplicon_min_depth=0, filter=None):
    """
        Extract major allele information from a JSNPmania variation line
        > _extract_major_snv_allele("sample1", ref, line2, nc_to_chr,1, major_allele_filter)
        >> ref = {'NC_000001.10': {}, 'NC_000017.10': {}, 'NC_000004.11': {}}
        >> line3 =  '4199\\t1711\\t0.40747797094546323\\tNC_000004.11\\t55968053\\tA\\t1711|0|2487|0\\tC\\tA\\t0.5924249642686994\\t428|426|426|431\\t0|0|0|0\\t654|590|575|668\\t0|0|0|0\\tchr4:55968023-55968093:+:384#chr4:55967974-55968185:+:243#chr4:55967975-55968186:-:188#chr4:55967939-55968092:+:232#chr4:55967938-55968091:-:314#chr4:55968027-55968093:-:333|0|chr4:55968023-55968093:+:549#chr4:55967974-55968185:+:321#chr4:55967975-55968186:-:307#chr4:55967939-55968092:+:437#chr4:55967938-55968091:-:413#chr4:55968027-55968093:-:436|0'
        >> line2 = '20\\t15\\t0.75\\tNC_000017.10\\t7579596\\tG\\t4|15|0|1\\tG\\tA\\t0.7894736842105263\\t2|0|0|2\\t14|0|0|1\\t0|0|0|0\\t1|0|0|0\\t0|chr17:7579592-7579839:+:13|0|0'
        >> line = '86\\t0\\t0.0\\tNC_000001.10\\t115252609\\tG\\t0|0|86|0\\tC\\tA/G/T\\t1.0\\t0|0|0|0\\t0|0|0|0\\t0|41|0|45\\t0|0|0|0\\t0|0|chr1:115252332-115252647:+:45#chr1:115252332-115252648:-:41|0'
        >> ref = _extract_ref_variant_info(ref,line3,1)
        >> ref = _extract_ref_variant_info(ref,line2,1)
        >> nc_to_chr = {'NC_000001.10': '1', 'NC_000017.10': '17', 'NC_000004.11': '4'}
        >> major_allele_filter = lambda allele_depth, columns: int(columns[0]) >= 20 and int(allele_depth)/int(columns[0]) >= 0.01 and (1.0 - float(columns[_column_converter_snv['allele_ratio']])) >= 0.01

        >> _extract_major_snv_allele("sample1", ref, line3, nc_to_chr,1, major_allele_filter)


    """
    columns = re.split("\t", line.rstrip())
    var_data = lambda field: columns[_column_converter_snv[field]]
    alleles = {}
    try:
        ref_data = ref[var_data('nc_number')][var_data('position')]
        variants = _ref_position.copy()
        del variants[var_data('reference')]
        allele_counts = _extract_allele_counts_snv(columns)
        nc_number = var_data('nc_number')
        position = var_data('position')
        #import itertools
        #var_iterator = iter(variants)
        #major_allele = [next(var_iterator)]
        #for v in var_iterator:
        #    if allele_counts[variants[v]] > allele_counts[variants[major_allele[0]]]:
        #         major_allele = [v]
        #    elif allele_counts[variants[v]] == allele_counts[variants[major_allele[0]]]:
        #         major_allele.append(v)
        for key, allele in variants.items():
            allele_depth = var_data('alleles_depth').split("|")[allele]
            major_allele_amplicon_information = var_data('amplicon_information').split("|")[allele]
            #print(str(allele) + " " + allele_depth + " " +  var_data('alleles_depth') + " " + str(variants.items()) + "\n")
            if filter is None or filter(allele_depth, columns):
                major_alle_info = "\t".join([
                    nc_to_chr[nc_number],
                    position,
                    position,
                    columns[_column_converter_snv['reference']],
                    key,
                    "comments: sample=" + sample +
                            " variantAlleleRatio=" + ("%.16f" %  (int(allele_depth)/(int(var_data('depth'))))).rstrip('0') +
                        " alleleFreq=" + str(ref_data['depth']) + "," + allele_depth +
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
                #alleles[nc_number + '#' + position + "#" + position + "#" + key] = major_alle_info
                alleles[nc_number + '#' + position + "#" + position + "#" + columns[_column_converter_snv['reference']] + "#" + key] = major_alle_info
    #except KeyError:
    #    logger.info("Unexpected reference base detected" + var_data('reference') + " when extracting major allele!")
    except AttributeError:
        logger.info("Unexpected reference base detected" + var_data('reference') + " when extracting major allele!")
    return alleles

def _extract_allele_counts_snv(columns):
    return [int(allele_d) for allele_d in columns[_column_converter_snv['alleles_depth']].split("|")]

def _extract_deletions(sample, ref, line, nc_to_chr, amplicon_min_depth=0, filter=None):
    """
        Extract information from a JSNPmania deletion line
    """
    columns = re.split("\t", line.rstrip())
    var_data = lambda field: columns[_column_converter_indel[field]]
    ref_data = ref[var_data('nc_number')][var_data('position')]
    position = int(var_data('position'))
    nc_number = var_data('nc_number')
    deletions = re.split("\|", var_data('indel'))
    extracted_deletions = {}
    for deletion in deletions:
        if deletion != '0':
            depth, variant = re.search("^(\d+)\(([0-9,-]+)\)$", deletion).groups()
            if filter is None or filter(depth,columns):
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
                if amplicon_min_depth:
                    for amp_info in var_data('amplicon_information').split("|"):
                        if amp_info.startswith("(" + variant + ")"):
                            amp_p, amp_m = _count_amplicons(amp_info,amplicon_min_depth)
                            deletion_info += " Tumor_var_plusAmplicons=" + str(amp_p) + \
                                              " Tumor_var_minusAmplicons=" + str(amp_m) + \
                                              " Tumor_ref_plusAmplicons=" + str(ref_data['amp+']) + \
                                              " Tumor_ref_minusAmplicons=" + str(ref_data['amp-']) + \
                                              " Tumor_var_ampliconInfo=" + amp_info + \
                                              " Tumor_ref_ampliconInfo=" + ref_data['ampInfo']
                            break;
                extracted_deletions[key] = deletion_info
    return extracted_deletions

def _extract_insertion(sample, ref, line, nc_to_chr, amplicon_min_depth=0, filter=None):
    """
        Extract insert information from a JSNPmania insertion line
    """
    columns = re.split("\t", line.rstrip())
    var_data = lambda field: columns[_column_converter_indel[field]]
    ref_data = ref[var_data('nc_number')][var_data('position')]
    insertions = re.split("\|", var_data('indel'))
    position = var_data('position')
    nc_number = var_data('nc_number')
    extracted_insertions = {}
    ref_depth = ref[nc_number][position]['depth']
    for index, insertion in enumerate(insertions):
        if insertion != '0':
            depth, variant = re.search("^(\d+)([A-Z]+$)", insertion).groups()
            if filter is None or filter(depth,columns):
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

def _update_deletion_data(variants, deletions):
    """
        >>> variants = {'NC_000007.13#140508151#140508154#-': 'ff', 'NC_000007.13#140508151#140508153#-': '7\\t140508151\\t140508153\\tCAA\\t-\\tcomments: sample=sample1 variantAlleleRatio=0.01646090534979424 alleleFreq=221,4 readDepth=243 Tumor_Del=0|15|0|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=(0,2):chr7:140508026-140508255:-:4 Tumor_ref_ampliconInfo=chr7:140508000-140508412:+:26#chr7:140508001-140508412:-:17#chr7:140508026-140508255:-:172'}
        >>> found_deletion = {'NC_000007.13#140508151#140508153#-': '7\\t140508151\\t140508153\\tCAG\\t-\\tcomments: sample=sample1 variantAlleleRatio=0.01656090534979424 alleleFreq=221,4 readDepth=243 Tumor_Del=0|15|0|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=(0,2):chr7:140508026-140508255:-:4 Tumor_ref_ampliconInfo=chr7:140508000-140508412:+:26#chr7:140508001-140508412:-:17#chr7:140508026-140508255:-:172'}
        >>> _update_deletion_data(variants, found_deletion)
        {'NC_000007.13#140508151#140508154#-': 'ff', 'NC_000007.13#140508151#140508153#-': '7\\t140508151\\t140508153\\tCAG\\t-\\tcomments: sample=sample1 variantAlleleRatio=0.01656090534979424 alleleFreq=221,4 readDepth=243 Tumor_Del=0|15|0|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=(0,2):chr7:140508026-140508255:-:4 Tumor_ref_ampliconInfo=chr7:140508000-140508412:+:26#chr7:140508001-140508412:-:17#chr7:140508026-140508255:-:172'}
    """
    for key, info in deletions.items():
        if key in variants:
            extract_info = lambda info: dict(tuple(comment.split("=")) for comment in info.split("\t")[-1].split(" ") if not comment.startswith("comment"))
            existing_variant_info = extract_info(variants[key])
            new_variant_info = extract_info(info)
            if float(new_variant_info['variantAlleleRatio']) > float(existing_variant_info['variantAlleleRatio']):
                variants[key] = info
        else:
            variants[key] = info
    return variants

def _get_amplicon_information(ref, nc_number, position, columns, ref_index, amplicon_min_depth):
    try:
        amplicon_information = columns[_column_converter_snv['amplicon_information']].split("|")[ref_index]
        amp_p, amp_m = _count_amplicons(amplicon_information, amplicon_min_depth)
        return (amplicon_information, amp_p, amp_m)
    except ValueError:
        return ('0', 0, 0)

def _extract_ref_variant_info(ref, line, amplicon_min_depth=0):
    """
        Extract reference variant information from a jSNPmania variations file.

        Parameters:
        ref: (dict) that will be populated with reference information, must already
            contain keys for the expected NC numbers, ex {'NC_000001.10': {}, 'NC_000002.8': {}}
        line: (string) entry from jSNPmania variations file
        amplicon_min_depth: (int) minum read depth for an amplicon if's going to be
            included in result, if set to null no amplicon information will be
            extracted.

        Returns: ref= {
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

        >>> ref = {'NC_000001.10': {}}
        >>> entry = '465\\t464\\t0.9978494623655914\\tNC_000001.10\\t162722741\\tG\\t1|464|0|0\\tG\\tA\\t0.9978494623655914\\t1|0|0|0\\t195|4|265|0\\t0|0|0|0\\t0|0|0|0\\tchr1:162722740-162722971:+:1|chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195|0|0'
        >>> _extract_ref_variant_info(ref, entry)
        {'NC_000001.10': {'162722741': {'reference': 'G', 'depth': 464}}}
        >>> _extract_ref_variant_info(ref, entry, 5)
        {'NC_000001.10': {'162722741': {'reference': 'G', 'depth': 464, 'amp+': 1, 'amp-': 1, 'ampInfo': 'chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195'}}}
        >>> _extract_ref_variant_info(ref, entry, 200)
        {'NC_000001.10': {'162722741': {'reference': 'G', 'depth': 464, 'amp+': 0, 'amp-': 1, 'ampInfo': 'chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195'}}}
        >>> entry = '8610\\t37\\t0.0042973286875725904\\tNC_000002.11\\t29416572\\tT\\t9|12|8552|37\\tC\\tT\\t0.9956921643963209\\t3|1|1|4\\t5|1|1|5\\t1863|2436|2415|1838\\t7|12|12|6\\tchr2:29416485-29416577:+:1|0|chr2:29416539-29416696:+:840#chr2:29416543-29416699:-:870#chr2:29416485-29416577:+:1496#chr2:29416485-29416577:-:2966#chr2:29416548-29416710:+:1322#chr2:29416548-29416710:-:979|0'
        >>> ref = {'NC_000002.11': {}}
        >>> _extract_ref_variant_info(ref, entry, 20)
        {'NC_000002.11': {'29416572': {'reference': 'T', 'depth': 37, 'amp+': 0, 'amp-': 0, 'ampInfo': '0'}}}
    """
    columns = re.split("\t", line.rstrip())
    try:
        nc_number = columns[_column_converter_snv['nc_number']]
        position = columns[_column_converter_snv['position']]
        ref_index = _ref_position[columns[_column_converter_snv['reference']]]
        ref[nc_number][position] = { 'reference': columns[_column_converter_snv['reference']],
                'depth': _extract_allele_counts_snv(columns)[ref_index]}
        if amplicon_min_depth:
            try:
                amp_info, amp_p, amp_m = _get_amplicon_information(ref,nc_number, position, columns,ref_index, amplicon_min_depth)
            except ValueError:
                amp_info, amp_p, amp_m = "0", 0, 0
            ref[nc_number][position]['amp+'] = amp_p
            ref[nc_number][position]['amp-'] = amp_m
            ref[nc_number][position]['ampInfo'] = amp_info
    except KeyError:
        logger.info("Unexpected reference base detected:" + columns[_column_converter_snv['reference']] + "!")
    return ref

def _count_amplicons(amplicon_information, min_depth=0):
    """
        Count number of HaloPlex amplicons for each strand that have fulfilled
        the minimum read depth.

        Parameters:
        amplicon_information: (string) a jSNPmania generate amplicon information
            string.
            SNV example
            chr1:162722740-162722971:+:1|chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195
            Insertion/Deletion example
            (0,0):chr1:162722740-162722971:+:1
        min_depth: (int) min required read depth for an amplicon to be included

        Returns: (int,int) --> (number of plus strand amplicons, number of minus strand amplicons)

        >>> _count_amplicons("chr1:162722740-162722971:+:1#chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195",2)
        (1, 1)
        >>> _count_amplicons("chr1:162722740-162722971:+:1#chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195")
        (2, 1)
        >>> _count_amplicons("(0,0):chr1:162722740-162722971:+:1")
        (1, 0)
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

if __name__ == "__main__":
    import doctest
    #import logging
    #import sys
    #logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    doctest.testmod()
    #result = doctest.testmod()
    #sys.exit(result.failed)
