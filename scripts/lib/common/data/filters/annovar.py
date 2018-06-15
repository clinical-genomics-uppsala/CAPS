def depth_and_vaf(depth, vaf, columns, headermap):
    return (int(columns[headermap['Read_depth']]) >= depth and
           float(columns[headermap['Variant_allele_ratio']]) >= vaf)

def overlap_region(columns, regions, headermap):
    entry_chr = columns[headermap['Chr']]
    entry_start = columns[headermap['Start']]
    entry_end = columns[headermap['End']]
    for (region_chr,region_start,region_end, *tail) in regions:
        if entry_chr == region_chr and entry_start <= region_end and entry_end >= region_start:
            return True
    else:
        return False

def is_exonic(columns, headermap):
    return columns[headermap['Type']].startswith('exonic')

def is_splicing(columns, headermap):
    return 'splicing' in columns[headermap['Type']]

def is_nonsynonymous(columns, headermap):
    return not columns[headermap['Exonic_type']].startswith("synonymous")

def in_1000g(vaf, columns, headermap):
    try:
        return vaf <= float(columns[headermap['Ratio_in_1000Genome']])
    except ValueError:
        return False

def check_key_exists(key, data):
    return key in data

def is_pindel_line(columns, headermap):
    return '+' in columns[headermap['Strands_Ins']] or '+' in columns[headermap['Strands_Del']]

def contains_valid_information(validate_levels, columns, headermap, ampliconmapped):
    if is_pindel_line(columns, headermap):
        return True
    else:
        try:
            var_plus = int(columns[headermap['#variant_+_amplicons']])
        except:
            var_plus = 0
        try:
            var_minus = int(columns[headermap['#variant_-_amplicons']])
        except:
            var_minus = 0
        try:
            ref_plus = int(columns[headermap['#reference_+_amplicons']])
        except:
            ref_plus = 0
        try:
            ref_minus = int(columns[headermap['#reference_-_amplicons']])
        except:
            ref_minus = 0
        try:
            ref = columns[headermap['Reference_base']]
        except:
            ref = "-"
        try:
            var = columns[headermap['Variant_base']]
        except:
            var = "-"
        return validate_levels(ref_plus, ref_minus, var_plus, var_minus, ref, var, ampliconmapped)
