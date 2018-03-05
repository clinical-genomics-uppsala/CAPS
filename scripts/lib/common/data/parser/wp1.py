
_column_converter_hotspot = {
    'chrom': 0,
    'start': 1,
    'end': 2,
    'gene': 3,
    'cds': 4,
    'aa': 5,
    'report': 6,
    'comment': 7,
    'exon': 8,
    'acc_num': 9,
    'amplicon_information': 10}

_valid_report_types = [
    'hotspot', # ToDo add description
    'indel', # ToDo add description
    'region', # ToDo add description
    'region_all', # ToDo add description
    'check' # ToDo add description
]

_column_converter_blacklist = {
    'chrom': 0,
    'start': 1,
    'stop': 2,
    'ref': 3,
    'var': 4,
    'gene': 5}

_column_converter_transcipt = {
    'gene': 0,
    'acc_num': 1,
    'protein_num': 2}

_column_converter_multibp_variants = {
    'chrom': 0,
    'start': 1,
    'stop': 2,
    'ref': 3,
    'var': 4,
    'gene': 5,
    'cds_change': 6,
    'aa_change': 7,
    'transcript': 8}

def valid_amplicon_nonhotspot(ref_plus, ref_minus, var_plus, var_minus):
    if (ref_plus + ref_minus) >= 3:
        if (var_plus + var_minus) >= 2:
            return True
    elif (ref_plus + ref_minus) >= 1:
        if (var_plus + var_minus) >= 1:
            return True
    elif (ref_plus + ref_minus) >= 0:
        if (var_plus + var_minus) >= 0:
            return True
    return False


def valid_amplicon_hotspot(ref_plus, ref_minus, var_plus, var_minus):
    if (var_plus == 0 and var_minus == 0):
        return False
    else:
        return True

class RejectingDict(dict):
    def __setitem__(self, k, v):
        if k in self.keys():
            raise ValueError("Key ({}) is already present!".format(k))
        else:
            return super(RejectingDict, self).__setitem__(k, v)

def remove_new_line(func):
    def func_wrapper(line):
        return func(line.rstrip('\n').rstrip('\r'))
    return func_wrapper

def skip_header(func):
    def func_wrapper(line):
        if line.startswith('#') or len(line) == 0:
            return None
        else:
            return func(line)
    return func_wrapper

def generate_key(report,chrom,start,stop):
    return (report,chrom,start,stop)


@skip_header
@remove_new_line
def parse_reference_info_line(line):
    columns = line.split("\t")
    return {columns[0].replace("chr",""): columns[1]}

def parse_reference_info(input_file):
    data = RejectingDict()
    def info_updater(key, value):
        data[key] = value
    with open(input_file, 'r') as info_lines:
        process_file(
            info_lines,
            parse_reference_info_line,
            [info_updater]
        )
    return data

@skip_header
@remove_new_line
def parse_hotspot_line(line):
    columns = line.split("\t")
    if not columns[_column_converter_hotspot['report']] in _valid_report_types:
        raise Exception("Unexpected report type {} provided!!!".format(columns[_column_converter_hotspot['report']]))
    start = columns[_column_converter_hotspot['start']]
    stop = columns[_column_converter_hotspot['end']]
    key = generate_key(columns[_column_converter_hotspot['report']],
                columns[_column_converter_hotspot['chrom']],
                start,
                stop)
    if columns[_column_converter_hotspot['comment']] == "-":
        columns[_column_converter_hotspot['comment']] = None
    return {key:
                {
                    'gene': columns[_column_converter_hotspot['gene']],
                    'cds': columns[_column_converter_hotspot['cds']],
                    'aa': columns[_column_converter_hotspot['aa']],
                    'comment': columns[_column_converter_hotspot['comment']],
                    'exon': columns[_column_converter_hotspot['exon']],
                    'acc_num': columns[_column_converter_hotspot['acc_num']],
                    'bwa': None,
                    'pindel': None,
                    'rd': ['NA'] * (int(stop) - int(start) + 1) if not columns[_column_converter_hotspot['report']] == 'indel' else 'NA',
                }
            }

def read_hotspot_file(input_file):  #ToDo make chromosom conversion
    hotspot_data = RejectingDict()
    intronic_data = RejectingDict()
    def hotspot_updater(key, value):
        if key[0] in hotspot_data:
                hotspot_data[key[0]][key[1:]] = value
        else:
            hotspot_data[key[0]] = {key[1:]: value}
    def indel_updater(key, value):
        if key[0] in hotspot_data and "intronic" in value['exon']:
            intronic_data[key[1:]] =value
    with open(input_file) as hotspot_lines:
        process_file(
            hotspot_lines,
            parse_hotspot_line,
            [hotspot_updater, indel_updater]
        )
        #for hotspotline in hotspot_lines:
        #    data = parse_hotspot_line(hotspot_line)
        #    if not data is None:
        #        for key, value in data.items():
        #            if key in hotspot:
        #                raise Exception("hotspot key {} already imported!!!".format(key))
        #            else:
        #                hotspot_data[key] = value
        #                if "intronic" in value['exon']:
        #                    if key in intronic_data:
        #                        raise Exception("intronic hotspot key {} already imported!!!".format(key))
        #                    else:
        #                        intronic_data[key] = value
    return (hotspot_data, intronic_data)

#def extract_hotspot_pindel_info(info, column_mapper, transcripts, multibp_data):
#    """
#        param: info: list with pindel variants
#    """
#    if info is None:
#        return info
#    for variant in info:

def get_read_level(read_levels, rd):
    try:
        for (level, depth_status, analyzable) in read_levels:
            rd = int(rd)
            if rd > int(level):
                return depth_status, analyzable
    except:
        pass
    return "-", "zero"

def get_transcript_info(variant, column_mapper, transcripts):
    aa, cds, acc_num, exon, exonic_type, commment = get_transcript_info(variant, column_mapper, transcripts)
    key = (
        info[column_mapper['Chr']],
        info[column_mapper['Start']],info[column_mapper['End']],
        info[column_mapper['Reference_base']],
        info[column_mapper['Variant_base']])
    # If variant have been defined in the multiple bp file, update aa, cds and acc_num with our predifiend values
    if key in multibp_data:
        aa = multibp_data[key]['aa']
        cds = multibp_data[key]['cds']
        if not acc_num == multibp_data[key]['nm']:
            comment = "altTranscript"
            acc_num = multibp_data[key]['nm']
    return aa, cds, acc_num, exon, exonic_type, comment

def merge_comment(comment1, comment2):
    if comment1 is None and comment2 is None:
        return "-"
    return " ".join([comment for comment in [comment1, comment2] if comment is not None ])

def bwa_and_pindel_overlap(var, pindel_variants, annovar_mapper):
    try:
        if pindel_variants['pindel'] is not None:
            for pindel in pindel_variants['pindel']:
                if int(pindel[annovar_mapper['Start']]) == int(var[annovar_mapper['Start']]) and \
                    int(pindel[annovar_mapper['End']]) == int(var[annovar_mapper['End']]):
                    return True
    except KeyError:
        pass
    return False;

def get_mutation_type(ref, var, insertion, deletion, mut_type):
    if ref == '-' or var == '-':
        return "2-indel"
    elif "+" in insertion or "+" in deletion:
        return "2-indel"
    return mut_type

def no_variant_found(sample, info, key, read_levels, annovar_mapper, mutation_type):
    variants =[]
    (entry_chr, entry_start, entry_end) = key
    for position in range(0, int(entry_end) - int(entry_start) + 1):
        found = level = read_depth = "NA"
        if isinstance(info['rd'], list):
            depth = info['rd'][position]
            try:
                if depth == "-":
                    found = "not in design"
                    level = "-"
                    read_depth = "-"
                else:
                    (level, found) = get_read_level(read_levels, depth)
                    if found == "yes":
                        found = "no"
                    read_depth = depth
            except KeyError as e:
                found = "not in design"
                level = "-"
                read_depth = "-"

        comment = info['comment']
        if comment is None:
            comment = "-"

        variants.append([
            sample,
            info['gene'],
            "-",
            info['exon'],
            info['aa'],
            info['cds'],
            info['acc_num'],
            comment,
            mutation_type,
            found,
            level,
            str(depth),
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            entry_chr,
            entry_start,
            entry_end,
            "-",
            "-",
            "-"])
    return variants

def extract_hotspot_snv_info(info, read_levels, annovar_mapper, transcripts, amplicon_mapped, mutation_type=None):
    try:
        if info['bwa'] is not None:
            variants = []
            for variant in info['bwa']:
                if not bwa_and_pindel_overlap(variant, info, annovar_mapper):
                    (aa, cds, acc_num, exon, exonic_type, comment) = get_transcript_info(variant, annovar_mapper, transcripts)
                    (found, level) = get_read_level(read_levels, info['read_depth'][0])
                    ref_plus, ref_minus, var_plus, var_minus, ref_all, var_all = get_amplicon_info(variant, annovar_mapper, amplicon_mapped)

                    # Make sure that amplicon mapped data contain the required number of amplicons
                    if amplicon_mapped and not valid_amplicon_hotspot(ref_plus, ref_minus, var_plus, var_minus):
                        found = "no"

                    # Use hotspot data if no transcript info can be extracted from annovar data
                    if exon == "-":
                        exon = info['exon']
                    if mutation_type is None:
                        mutation_type = get_mutation_type(variant[annovar_mapper['Reference_base']],variant[annovar_mapper['Variant_base']],
                                        variant[annovar_mapper['Strands_Ins']],variant[annovar_mapper['Strands_Del']],mutation_type)

                    comment = merge_comment(info['comment'],comment)

                    variants.append(
                        [
                            variant[annovar_mapper['Sample']],
                            variant[annovar_mapper['Gene']],
                            exonic_type,
                            exon,
                            aa,
                            cds,
                            acc_num,
                            comment,
                            mutation_type,
                            level,
                            found,
                            variant[annovar_mapper['Read_depth']],
                            variant[annovar_mapper['#reference_alleles']],
                            variant[annovar_mapper['#_variant_alleles']],
                            variant[annovar_mapper['Variant_allele_ratio']],
                            variant[annovar_mapper['dbSNP_id']],
                            variant[annovar_mapper['Ratio_in_1000Genome']],
                            variant[annovar_mapper['ESP_6500']],
                            variant[annovar_mapper['Clinically_flagged_dbSNP']],
                            variant[annovar_mapper['Cosmic']],
                            variant[annovar_mapper['ClinVar_CLNDBN']],
                            variant[annovar_mapper['ClinVar_CLINSIG']],
                            ref_plus,
                            ref_minus,
                            var_plus,
                            var_minus,
                            variant[annovar_mapper['Strands_A']],
                            variant[annovar_mapper['Strands_G']],
                            variant[annovar_mapper['Strands_C']],
                            variant[annovar_mapper['Strands_T']],
                            variant[annovar_mapper['Strands_Ins']],
                            variant[annovar_mapper['Strands_Del']],
                            ref_all,
                            var_all,
                            variant[annovar_mapper['Chr']],
                            variant[annovar_mapper['Start']],
                            variant[annovar_mapper['End']],
                            variant[annovar_mapper['Reference_base']],
                            variant[annovar_mapper['Variant_base']],
                            variant[annovar_mapper['Transcripts']]
                        ]
                    )
            return variants
    except KeyError as e:
        if str(e) != "bwa":
            raise e
    return None#print("Error " + str(e))
        #return None

def extract_hotspot_pindel_info(info, read_levels, annovar_mapper, transcripts, amplicon_mapped, mutation_type=None):
    try:
        if info['pindel'] is not None:
            variants = []
            for variant in info['pindel']:
                (aa, cds, acc_num, exon, exonic_type, comment) = get_transcript_info(variant, annovar_mapper, transcripts)
                (found, level) = get_read_level(read_levels, info['read_depth'][0])
                ref_plus, ref_minus, var_plus, var_minus, ref_all, var_all = get_amplicon_info(variant,annovar_mapper, amplicon_mapped)

                # Make sure that amplicon mapped data contain the required number of amplicons
                if amplicon_mapped and not valid_amplicon_hotspot(ref_plus, ref_minus, var_plus, var_minus):
                    found = "no"
                # Use hotspot data if no transcript info can be extracted from annovar data
                if exon == "-":
                    exon = info['exon']
                comment = merge_comment(info['comment'],comment)

                if mutation_type is None:
                    mutation_type = get_mutation_type(variant[annovar_mapper['Reference_base']],variant[annovar_mapper['Variant_base']],
                                    variant[annovar_mapper['Strands_Ins']],variant[annovar_mapper['Strands_Del']],mutation_type)
                variants.append(
                    [
                        variant[annovar_mapper['Sample']],
                        variant[annovar_mapper['Gene']],
                        exonic_type,
                        exon,
                        aa,
                        cds,
                        acc_num,
                        comment,
                        mutation_type,
                        level,
                        found,
                        variant[annovar_mapper['Read_depth']],
                        variant[annovar_mapper['#reference_alleles']],
                        variant[annovar_mapper['#_variant_alleles']],
                        variant[annovar_mapper['Variant_allele_ratio']],
                        variant[annovar_mapper['dbSNP_id']],
                        variant[annovar_mapper['Ratio_in_1000Genome']],
                        variant[annovar_mapper['ESP_6500']],
                        variant[annovar_mapper['Clinically_flagged_dbSNP']],
                        variant[annovar_mapper['Cosmic']],
                        variant[annovar_mapper['ClinVar_CLNDBN']],
                        variant[annovar_mapper['ClinVar_CLINSIG']],
                        ref_plus,
                        ref_minus,
                        var_plus,
                        var_minus,
                        variant[annovar_mapper['Strands_A']],
                        variant[annovar_mapper['Strands_G']],
                        variant[annovar_mapper['Strands_C']],
                        variant[annovar_mapper['Strands_T']],
                        variant[annovar_mapper['Strands_Ins']],
                        variant[annovar_mapper['Strands_Del']],
                        ref_all,
                        var_all,
                        variant[annovar_mapper['Chr']],
                        variant[annovar_mapper['Start']],
                        variant[annovar_mapper['End']],
                        variant[annovar_mapper['Reference_base']],
                        variant[annovar_mapper['Variant_base']],
                        variant[annovar_mapper['Transcripts']]
                    ]
                )
            return variants
    except KeyError as e:
        if str(e) != "pindel":
            raise e
    #    else:
    #        print("Error " + str(e))
    return None

@skip_header
@remove_new_line
def parse_blacklist_line(line):
    columns = line.split("\t")
    key = ('blacklist',
                columns[_column_converter_blacklist['chrom']],
                columns[_column_converter_blacklist['start']],
                columns[_column_converter_blacklist['stop']],
                columns[_column_converter_blacklist['ref']],
                columns[_column_converter_blacklist['var']])
    return {key: {'gene': columns[_column_converter_blacklist['gene']]}}

def read_blacklist(input_file): #ToDo make chromosom conversion
    blacklist_data = RejectingDict()
    def blaclist_updater(key,value):
        blacklist_data[key[1:]] = value
    with open(input_file) as blacklist_lines:
        process_file(
            blacklist_lines,
            parse_blacklist_line,
            [blaclist_updater])
        #for blacklist_line in blacklist_lines:
        #    data = parse_blacklist_line(blacklist_line)
        #    if not data is None:
        #        for key, value in data.items():
        #            if key in blacklist_data:
        #                raise Exception("blacklist key {} already imported!!!".format(key))
        #            else:
        #                blacklist_data[key] = value
    return blacklist_data

@skip_header
@remove_new_line
def parse_transcript_line(line):
    columns = line.split("\t")
    key = ('maintranscript',
            columns[_column_converter_transcipt['gene']],
            columns[_column_converter_transcipt['acc_num']].split(".")[0])
    return {key: {'gene': columns[_column_converter_transcipt['gene']]}}  #ToDo maybe remove version

def read_maintranscripts(input_file):
    maintranscript_data = RejectingDict()
    def transcript_updater(key,value):
        maintranscript_data[key[1]] = key[2]
    with open(input_file) as maintranscript_lines:
        process_file(
            maintranscript_lines,
            parse_transcript_line,
            [transcript_updater])
        #for maintranscript_line in maintranscript_lines:
        #    data = parse_maintranscript_line(maintranscript_line)
        #    if not data is None:
        #        for key, value in data.items():
        #            if key in maintranscript_data:
        #                raise Exception("maintranscript key {} already imported!!!".format(key))
        #            else:
        #                maintranscript_data[key] = value
    return maintranscript_data

@skip_header
@remove_new_line
def parse_multi_bp_variant_line(line):
    columns = line.split("\t")
    key = ('multibp',
            columns[_column_converter_multibp_variants['chrom']],
            columns[_column_converter_multibp_variants['start']],
            columns[_column_converter_multibp_variants['stop']],
            columns[_column_converter_multibp_variants['ref']],
            columns[_column_converter_multibp_variants['var']])
    return {key: {
                'cds': columns[_column_converter_multibp_variants['cds_change']],
                'aa': columns[_column_converter_multibp_variants['aa_change']],
                'nm': columns[_column_converter_multibp_variants['transcript']]}} #ToDo maybe remove version

def read_multibpvariants(input_file):
    multibp_variant_data = RejectingDict()#{}
    with open(input_file) as multibp_variant_lines:
        #for multibp_variant_line in multibp_variant_lines:
        #    data = parse_multi_bp_variant_line(multibp_variant_line)
        #    if not data is None:
        #        for key, value in data.items():
        #            if key in multibp_variant_data:
        #                raise Exception("multibp variant key {} already imported!!!".format(key))
        #            else:
        #                multibp_variant_data[key] = value
        def multibp_updater(key,value):
            multibp_variant_data[key] = value
        process_file(
            multibp_variant_lines,
            parse_multi_bp_variant_line,
            [multibp_updater])
    return multibp_variant_data

def process_file(file_handler, line_parser, data_updater):
    for line in file_handler:
        data = line_parser(line)
        if data is not None:
            for key, value in data.items():
                for updater in data_updater:
                    updater(key,value)

def get_amplicon_info(columns, column_mapper, ampliconMapped):
    if ampliconMapped:
        return (columns[column_mapper['#reference_+_amplicons']],
                columns[column_mapper['#reference_-_amplicons']],
                columns[column_mapper['#variant_+_amplicons']],
                columns[column_mapper['#variant_-_amplicons']],
                columns[column_mapper['Reference_ampliconinfo']],
                columns[column_mapper['Variant_ampliconinfo']])
    return ('NA','NA','NA','NA','NA','NA')

def get_transcript_info(columns, column_mapper, transcripts):
    exonic_type = columns[column_mapper['Exonic_type']]
    gene = columns[column_mapper['Gene']]
    if "splicing" in columns[column_mapper['Type']]:
        exonic_type = columns[column_mapper['Type']]
    aa = cds = acc_num = exon = "-"  # set the parameters to - as default
    comm = None  # If variable to add comment in
    found = False;
    if not "-" == columns[column_mapper['Transcripts']]:  # Check that there are any transcript
        allTranscript = columns[column_mapper['Transcripts']].split(",")
        # If the gene exist in preferd transcript hash use that transcript, otherwise take the first
        if gene in transcripts:  # check if the genename exist in the transcript hash
            for tr in allTranscript:  # if so go through and see if any transcript overlaps with the main transcript in hash
                if transcripts[gene] in tr:
                    transcriptInfo = tr.split(":")
                    found = True  # The main transript is found in the list of transcript for this mutation
                    if len(transcriptInfo) >= 5:
                        aa = transcriptInfo[4]
                        cds = transcriptInfo[3]
                        exon = transcriptInfo[2]
                        acc_num = transcriptInfo[1]
                    elif len(transcriptInfo) >= 4:
                        cds = transcriptInfo[3]
                        exon = transcriptInfo[2]
                        acc_num = transcriptInfo[1]
                    elif len(transcriptInfo) >= 3:
                        exon = transcriptInfo[2]
                        acc_num = transcriptInfo[1]
                    elif len(transcriptInfo) >= 2:
                        acc_num = transcriptInfo[1]
            if not found:
                comm = "altTranscript"
        else:
            transcriptInfo = allTranscript[0].split(":")

            if len(transcriptInfo) >= 5:
                aa = transcriptInfo[4]
                cds = transcriptInfo[3]
                exon = transcriptInfo[2]
                acc_num = transcriptInfo[1]
            elif len(transcriptInfo) >= 4:
                cds = transcriptInfo[3]
                exon = transcriptInfo[2]
                acc_num = transcriptInfo[1]
            elif len(transcriptInfo) >= 3:
                exon = transcriptInfo[2]
                acc_num = transcriptInfo[1]
            elif len(transcriptInfo) >= 2:
                acc_num = transcriptInfo[1]

    return aa, cds, acc_num, exon, exonic_type, comm
