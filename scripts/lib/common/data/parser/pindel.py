# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#def convert_to_annovar_input(sample_name, output, deletions_file, insertions_file, nc_to_chr, min_allele_ratio, min_read_depth):
#    with open(deletions) as deletions:

import vcf

_column_converter_pindel = {'indel_info': 2 ,'chr': 3, 'start': 4, 'stop': 5, 'read_info': 19}
_pindel_read_info_column = {'sample_name': 0,
                            'left_breakpoint_read_depth': 1,
                            'right_breakpoint_read_depth': 2,
                            'reads_alternate_allele_fortward': 3,
                            'uniq_reads_alternate_allele_fortward': 4,
                            'reads_alternate_allele_reverse': 5,
                            'uniq_reads_alternate_allele_reverse': 6}
def __extract_value(info):
    return info.split()[-1]

def _extract_deletion(line, amplicon_min_depth=0, filter=None):
    columns = re.split(r'\s', line.rstrip())
    pindel_data = lambda field: columns[_column_converter_indel[field]]
    chr = pindel_data('chr').split(r'\s')[1]
    start = int(pindel_data('start').split(r'\s')[1])
    end = start - 1
    insert_info = pindel_data('insert_info').split(r'\s')

def parse_pindel_line(line,report_min_depth=False):
    """
        Function used to parse lines from pindel output.

    >>> line_del = \
            '0\\tD 4\\tNT 0 ""\\tChrID chr1\\tBP 67443969\\t67443974\\t' + \
            'BP_range 67443969\\t67443980\\tSupports 5\\t1\\t+ 5\\t1\\t- 0\\t0\\t' + \
            'S1 6\\tSUM_MS 300\\t1\\tNumSupSamples 1\\t1\\tSample_name 1 0 5 1 0 0'
    >>> line_ins = '0\\tI 51\\t' + \
            'NT 51 "CATATATGGGGCCAAAGGAACAACTCCATGTTTTCTCAAAGGCCTAGAGAA"\\t' + \
    	    'ChrID chr1\\tBP 31581378\\t31581379\\t' + \
            'BP_range 31581378\\t31581437\\tSupports 3\\t1\\t+ 0\\t0\\t- 3\\t1\\t' + \
            'S1 4\\tSUM_MS 159\\t1\\tNumSupSamples 1\\t1\\tSample_name 5 5 0 0 3 1'
    >>> parse_pindel_line(line_del)
    {'chr': 'chr1', 'start': 67443969, 'stop': 67443974, 'indel': '-', 'depth': 1, 'fwd': 5, 'fwd_u': 1, 'rve': 0, 'rve_u': 0}
    >>> parse_pindel_line(line_del, True)
    {'chr': 'chr1', 'start': 67443969, 'stop': 67443974, 'indel': '-', 'depth': 0, 'fwd': 5, 'fwd_u': 1, 'rve': 0, 'rve_u': 0}
    >>> parse_pindel_line(line_ins)
    {'chr': 'chr1', 'start': 31581378, 'stop': 31581379, 'indel': 'CATATATGGGGCCAAAGGAACAACTCCATGTTTTCTCAAAGGCCTAGAGAA', 'depth': 5, 'fwd': 0, 'fwd_u': 0, 'rve': 3, 'rve_u': 1}
    >>> parse_pindel_line(line_ins, True)
    {'chr': 'chr1', 'start': 31581378, 'stop': 31581379, 'indel': 'CATATATGGGGCCAAAGGAACAACTCCATGTTTTCTCAAAGGCCTAGAGAA', 'depth': 5, 'fwd': 0, 'fwd_u': 0, 'rve': 3, 'rve_u': 1}
    """
    columns = line.rstrip('\r\n').split('\t')
    # column = 'ChrID chr1'
    chr = __extract_value(columns[_column_converter_pindel['chr']])
    # column = 'BP_range 67443969'
    start = int(__extract_value(columns[_column_converter_pindel['start']]))
    # column = '31581379'
    stop = int(__extract_value(columns[_column_converter_pindel['stop']]))
    #column = 'NT 0 ""'
    #column = 'NT 5 "CATAT"'
    indel = columns[_column_converter_pindel['indel_info']].split()[2] \
                    .replace('""', "-").replace('"', "")

    # column = 'Sample_name 5 5 0 0 3 1'
    read_info = columns[_column_converter_pindel['read_info']].split()
    left_breakpoint_read_depth = \
        int(read_info[_pindel_read_info_column['left_breakpoint_read_depth']])
    right_breakpoint_read_depth = \
        int(read_info[_pindel_read_info_column['right_breakpoint_read_depth']])
    allternate_allele_reads_fwd = \
        int(read_info[_pindel_read_info_column['reads_alternate_allele_fortward']])
    allternate_allele_reads_rve = \
        int(read_info[_pindel_read_info_column['reads_alternate_allele_reverse']])
    allternate_allele_uniq_reads_fwd = \
        int(read_info[_pindel_read_info_column['uniq_reads_alternate_allele_fortward']])
    allternate_allele_uniq_reads_rve = \
        int(read_info[_pindel_read_info_column['uniq_reads_alternate_allele_reverse']])
    read_depth = left_breakpoint_read_depth

    if report_min_depth and \
            left_breakpoint_read_depth >  right_breakpoint_read_depth:
        read_depth = right_breakpoint_read_depth
    elif not report_min_depth and \
            left_breakpoint_read_depth <  right_breakpoint_read_depth:
        read_depth = right_breakpoint_read_depth

    return {'chr': chr, 'start': start, 'stop': stop, 'indel': indel,
                'depth': read_depth,
                'fwd': allternate_allele_reads_fwd,
                'fwd_u': allternate_allele_uniq_reads_fwd,
                'rve': allternate_allele_reads_rve,
                'rve_u': allternate_allele_uniq_reads_rve}

def _get_start_pos(pos,type):
    if type == "INS":
        return pos
    elif type in ["DEL", "RPL"]:
        return pos + 1
    else:
        raise Exception("Unhandled indel: " + str(type))

def _extract_alt_allele_from_vcf(alt,type):
    if type == "SNV":
        return str(alt)
    elif type == "DEL":
        return "-"
    elif type in ["RPL", "INS"]:
        return str(alt)[1:]
    else:
        raise Exception("Unhandled indel: " + str(type) + " " + str(alt))

def _extract_ref_allele_from_vcf(alt,type):
    if type == "SNV":
        return str(alt)
    elif type == "INS":
        return "-"
    elif type in ["RPL", "DEL"]:
        return str(alt)[1:]
    else:
        raise Exception("Unhandled indel: " + str(type) + " " + str(alt))

# ToDO need to rename read_method
def convert_to_annovar_input(sample, output_file, input_vcf_file, pindel_deletions_file, pindel_insertions_file, min_read_depth = 1, min_vaf = 0, read_method = "min"):
    if not read_method in ["min", "max"]:
        raise Exception(read_method, "Expeced input \"min\" or \"max\"")
    def create_key(chr,start,stop,var): return chr + "#" + str(start) + "#" + str(stop) + "#" + var
    deletions = {}
    with open(output_file, 'w') as output:
        with open(pindel_deletions_file) as deletion_lines:
            for line in deletion_lines:
                if "Chr" in line:
                    data = parse_pindel_line(line,read_method == "min")
                    deletions[create_key(data['chr'],data['start'],data['stop'], data["indel"])] = data
        insertions = {}
        with open(pindel_insertions_file) as insertion_lines:
            for line in insertion_lines:
                if "Chr" in line:
                    data = parse_pindel_line(line, read_method == "min")
                    insertions[create_key(data['chr'],data['start'],data['stop'], data["indel"])] = data
        with open(input_vcf_file, 'r') as vcf_file:
            vcf_reader =  vcf.Reader(vcf_file)
            for record in vcf_reader:
                for alt in record.ALT:
                    key = create_key(record.CHROM, record.POS, record.INFO['END'] + 1, _extract_alt_allele_from_vcf(alt,record.INFO['SVTYPE']))
                    var_reads, ref_reads = 0, 0
                    strand_info = ""
                    if record.INFO['SVTYPE'] in ["DEL", "RPL"]:
                        if key in deletions:
                                strand_info = "Tumor_Del=+" + str(deletions[key]['fwd']) + "|-" + str(deletions[key]['rve'])
                                var_reads = int(deletions[key]['fwd']) +  int(deletions[key]['rve'])
                                ref_reads = int(deletions[key]['depth'])
                        else:
                            continue
                    elif record.INFO["SVTYPE"] == "INS":
                        if key in insertions:
                            strand_info = "Tumor_Ins=+" + str(insertions[key]['fwd']) + "|-" + str(insertions[key]['rve'])
                            var_reads = int(insertions[key]['fwd']) +  int(insertions[key]['rve'])
                            ref_reads = int(insertions[key]['depth'])
                        else:
                            continue
                    else:
                        print("Unhandled type: "  + str(record.INFO["SVTYPE"]) + " " +   + str(record))
                        continue
                    read_depth = var_reads + ref_reads
                    if var_reads / read_depth >= min_vaf and read_depth >= min_read_depth:
                        output.write(_generate_output(record, sample, alt, ref_reads, var_reads, strand_info) + "\n")


def _generate_output(record, sample, alt, ref_reads, var_reads, strand_info):
    return "\t".join([
        record.CHROM.replace("chr",""),
        str(_get_start_pos(record.POS,record.INFO['SVTYPE'])),
        str(record.INFO['END']),
        _extract_ref_allele_from_vcf(record.REF,record.INFO['SVTYPE']),
        _extract_alt_allele_from_vcf(alt,record.INFO['SVTYPE']),
        " ".join(["comments: sample=" + sample,
            "variantAlleleRatio=" + str(var_reads / (var_reads + ref_reads)),
            "alleleFreq=" + str(ref_reads) + "," + str(var_reads),
            "readDepth=" + str(ref_reads + var_reads),
            strand_info,
            "Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=- Tumor_ref_minusAmplicons=-"])])

if __name__ == "__main__":
    import doctest
    #import logging
    #import sys
    #logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    doctest.testmod()
    #result = doctest.testmod()
    #sys.exit(result.failed)
