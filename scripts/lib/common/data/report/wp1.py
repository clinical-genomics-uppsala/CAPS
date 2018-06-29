from scripts.lib.common.data.parser.wp1 import read_hotspot_file, read_blacklist, read_maintranscripts, read_multibpvariants, parse_reference_info, get_transcript_info,get_amplicon_info, extract_hotspot_pindel_info, extract_hotspot_snv_info
from scripts.lib.common.data.parser.wp1 import valid_indel, valid_amplicon_nonhotspot, valid_amplicon_hotspot, get_read_level, merge_comment, no_variant_found
from scripts.lib.common.data.parser.jsnpmania import _column_converter_snv, extract_ref_variant_info, _column_converter_snv
from scripts.lib.common.data.parser.annovar import generate_headermap
from scripts.lib.common.data.filters.annovar import depth_and_vaf, overlap_region, is_exonic
from scripts.lib.common.data.filters.annovar import is_splicing, is_nonsynonymous, in_1000g
from scripts.lib.common.data.filters.annovar import check_key_exists,contains_valid_information
from scripts.lib.common.data.filters.annovar import is_pindel_line
from scripts.lib.common.data.filters.wp1 import and_condition, or_condition, evalute_else, gene_name, insertion_var, vaf_above_or_equal, base_in_ref_or_var, evaluate_tuple
from scripts.lib.common.data.parser.common import create_chr_mapper, import_ampregion_seq, import_bed_file, region_data_structure_generator, _column_bed_file

def create_filtered_mutations_header():
    return "\t".join([
        "Sample",
        "Gene",
        "Variant_type",
        "Exon",
        "AA_change",
        "CDS_change",
        "Accession_number",
        "Comment",
        "Report",
        "Found",
        "Min_read_depth300",
        "Total_read_depth",
        "Reference_read_depth",
        "Variant_read_depth",
        "Variant_allele_ratio",
        "dbSNP_id",
        "Ratio_1000G",
        "Ratio_ESP6500",
        "Clinically_flagged_dbSNP",
        "Cosmic",
        "ClinVar_CLNDB",
        "Clinval_CLINSIG",
        "Reference_plus_amplicons",
        "Reference_minus_amplicons",
        "Variant_plus_amplicons",
        "Variant_minus_amplicons",
        "Strands_A_F+F-S+S-",
        "Strands_G_F+F-S+S-",
        "Strands_C_F+F-S+S-",
        "Strands_T_F+F-S+S-",
        "Strands_Ins",
        "Strands_Del",
        "Ref_aligned_amplicons",
        "Var_aligned_amplicons",
        "Chr",
        "Start",
        "End",
        "Reference_base",
        "Variant_base",
        "All_transcripts_annotation"
    ])

def filter_annovar_variant(intronic_data, columns, headermap, blacklist_data, chr, start, stop, key, min_read_depth, min_vaf,genome1000_vaf):
    # Filter data. Ignores variants below a certain read depth, vaf or if it's found in the 1000 genome
    if depth_and_vaf(min_read_depth, min_vaf, columns, headermap) and \
        not in_1000g(genome1000_vaf, columns, headermap):
        # Save data if it's intronic hotspot, else it must be synonymous
        # and exonic or splicing.
        if overlap_region(columns,intronic_data,headermap) or \
            (is_nonsynonymous(columns, headermap) and \
            (is_exonic(columns, headermap) or is_splicing(columns, headermap))):
            # Exclude, if it's found in the blacklist
            if not check_key_exists(key, blacklist_data):
                return True
    return False

def extract_msi_markers(input_file, output_file, tgfbr2_1bp = 0.01,tgfbr2_2bp = 0.01, acvr2a_1bp = 0.01, acvr2a_2bp = 0.01):
    def filter_rows(columns,mapper):
        if gene_name("TGFBR2")(columns,mapper):
            if insertion_var(1)(columns,mapper):
                return vaf_above_or_equal(tgfbr2_1bp)(columns,mapper)
            else:
                return vaf_above_or_equal(tgfbr2_2bp)(columns,mapper)
        elif gene_name("ACVR2A")(columns,mapper):
            if insertion_var(1)(columns,mapper):
                return vaf_above_or_equal(acvr2a_1bp)(columns,mapper)
            else:
                return vaf_above_or_equal(acvr2a_2bp)(columns,mapper)
        return False

    with open(output_file,'w') as output:
        output.write("#" + create_filtered_mutations_header())
        with open(input_file, 'r') as mutations:
            mapper = generate_headermap(mutations.readline(),"#Sample")
            for line in mutations:
                columns = line.rstrip("\n").rstrip("\r").split("\t")
                if filter_rows(columns, mapper):
                    output.write("\n" + "\t".join(columns))


def generate_filtered_mutations(sample, output_file, hotspot_file, snpmania_variants_file, annovar_files, multibp_file, chr_to_nc_file, min_read_depth, min_vaf, read_depth_classes, genome1000_vaf, ampliconmapped, blacklist_file, transcript_file):
    """
        >>> hotspot = '/snakemake/resources/Mutations_Lung_20171219.csv'
        >>> snpmania = '/snakemake/resources/18-215.ampliconmapped.variations'
        >>> annovar = '[/snakemake/resources/annovarOutput',/snakemake/resources/annovarOutput']
        >>> multibp = '/snakemake/resources/MultipleBp_variations.csv'
        >>> min_read_depth = 30
        >>> min_vaf = 0.01
        >>> read_depth_classes = [(300 ,"ok","yes"), (:q ,"low","yes"),(0 ,"low","not analyzable")]
        >>> max_1000genome = 0.02
        >>> blacklist = '/snakemake/resources/blacklist_variantClusterAllvarPindel_20151110.txt'
        >>> transcript = '/snakemake/resources/mainTranscripts.txt'
        >>> chr_to_nc = '/snakemake/resources/reference.chr_to_nc.hg19.info'
        >>> generate_filtered_mutations("sample1","/snakemake/resources/wp1.report", hotspot, snpmania, annovar, multibp, chr_to_nc, min_read_depth, min_vaf, read_depth_classes, max_1000genome, True, blacklist, transcript)
        "F"
    """
    chr_to_nc = parse_reference_info(chr_to_nc_file)
    (hotspot_data, intronic_data) = read_hotspot_file(hotspot_file)
    blacklist_data = read_blacklist(blacklist_file)
    transcript_data = read_maintranscripts(transcript_file)
    multibp_data = read_multibpvariants(multibp_file)
    annovar_header_mapper = None
    with open(output_file, 'w') as output:
        output.write("#" + create_filtered_mutations_header())
        for annovar_file in annovar_files:
            with open(annovar_file,'r') as annovar_lines:
                annovar_header_mapper = generate_headermap(annovar_lines.readline(),"#Sample")
                for line in annovar_lines:
                    columns = line.rstrip("\n").rstrip("\r").split("\t")
                    columns[annovar_header_mapper['Chr']] = chr_to_nc[columns[annovar_header_mapper['Chr']]]
                    key = (chrom, start, stop) = (columns[annovar_header_mapper['Chr']], columns[annovar_header_mapper['Start']], columns[annovar_header_mapper['End']])
                    blacklist_key = (columns[annovar_header_mapper['Chr']], columns[annovar_header_mapper['Start']], columns[annovar_header_mapper['End']], columns[annovar_header_mapper['Reference_base']], columns[annovar_header_mapper['Variant_base']])

                    if filter_annovar_variant(intronic_data, columns, annovar_header_mapper, blacklist_data, chr, start, stop, blacklist_key, min_read_depth, min_vaf,genome1000_vaf):
                        added = False
                        data_source = 'bwa'
                        #Should be able to remov check for source!
                        if is_pindel_line(columns, annovar_header_mapper):
                            data_source = 'pindel'
                        for (amplicon_validator, reports) in [(valid_amplicon_hotspot, ('hotspot',)), (valid_amplicon_nonhotspot ,('region','region_all')), (valid_indel,('indel',))]:
                            if not added:
                                if data_source == 'pindel' or contains_valid_information(amplicon_validator,columns, annovar_header_mapper, ampliconmapped):
                                    for report in reports:
                                        if not added:
                                            try:
                                                if key in hotspot_data[report]:
                                                    try:
                                                        hotspot_data[report][key][data_source].append(columns)
                                                    except (AttributeError, KeyError) as e:
                                                        hotspot_data[report][key][data_source] = [columns]
                                                    added = True
                                                    break
                                                else:
                                                    (chrom, start, end) = (columns[annovar_header_mapper['Chr']], columns[annovar_header_mapper['Start']], columns[annovar_header_mapper['End']])
                                                    for entry_key in hotspot_data[report]:
                                                        (entry_chrom, entry_start, entry_end) = entry_key
                                                        if entry_chrom == chrom and \
                                                            start <= entry_end and end >= entry_start:
                                                            try:
                                                                hotspot_data[report][entry_key][data_source].append(columns)
                                                            except (AttributeError, KeyError) as e:
                                                                hotspot_data[report][entry_key][data_source] = [columns]
                                                            added = True
                                                            break
                                            except KeyError:
                                                continue
                        if not added:
                            (ref_plus, ref_minus, var_plus, var_minus, ref_all, var_all) = get_amplicon_info(columns, annovar_header_mapper, ampliconmapped)

                            if data_source == 'pindel' or \
                                (contains_valid_information(valid_amplicon_nonhotspot, columns, annovar_header_mapper, ampliconmapped) or contains_valid_information(valid_indel, columns, annovar_header_mapper, ampliconmapped)):
                                (found, level) = get_read_level(read_depth_classes, columns[annovar_header_mapper['Read_depth']])
                                (aa, cds, acc_num, exon, exonic_type, comment) =  get_transcript_info(columns, annovar_header_mapper, transcript_data)
                                comment = merge_comment(comment,"-")
                                output.write("\n" + "\t".join([
                                    columns[annovar_header_mapper['Sample']],
                                    columns[annovar_header_mapper['Gene']],
                                    exonic_type,
                                    exon,
                                    aa,
                                    cds,
                                    acc_num,
                                    comment,
                                    "4-other",
                                    level,
                                    found,
                                    columns[annovar_header_mapper['Read_depth']],
                                    columns[annovar_header_mapper['#reference_alleles']],
                                    columns[annovar_header_mapper['#_variant_alleles']],
                                    columns[annovar_header_mapper['Variant_allele_ratio']],
                                    columns[annovar_header_mapper['dbSNP_id']],
                                    columns[annovar_header_mapper['Ratio_in_1000Genome']],
                                    columns[annovar_header_mapper['ESP_6500']],
                                    columns[annovar_header_mapper['Clinically_flagged_dbSNP']],
                                    columns[annovar_header_mapper['Cosmic']],
                                    columns[annovar_header_mapper['ClinVar_CLNDBN']],
                                    columns[annovar_header_mapper['ClinVar_CLINSIG']],
                                    ref_plus,
                                    ref_minus,
                                    var_plus,
                                    var_minus,
                                    columns[annovar_header_mapper['Strands_A']],
                                    columns[annovar_header_mapper['Strands_G']],
                                    columns[annovar_header_mapper['Strands_C']],
                                    columns[annovar_header_mapper['Strands_T']],
                                    columns[annovar_header_mapper['Strands_Ins']],
                                    columns[annovar_header_mapper['Strands_Del']],
                                    ref_all,
                                    var_all,
                                    columns[annovar_header_mapper['Chr']],
                                    columns[annovar_header_mapper['Start']],
                                    columns[annovar_header_mapper['End']],
                                    columns[annovar_header_mapper['Reference_base']],
                                    columns[annovar_header_mapper['Variant_base']],
                                    columns[annovar_header_mapper['Transcripts']]
                                    ]))
        intronic_data = {}
        with open(snpmania_variants_file, 'r') as variants:
            for line in variants:
                if not line.startswith("#"):
                    columns = line.rstrip("\n").rstrip("\r").split("\t")
                    nc_chr, position = columns[_column_converter_snv['nc_number']], columns[_column_converter_snv['position']]
                    for report in ['indel', 'region', 'region_all', 'hotspot']:
                        if report in hotspot_data:
                            for (report_nc, report_start, report_end ) in hotspot_data[report]:
                                if nc_chr == report_nc and int(report_start) <= int(position) <= int(report_end):
                                    try:
                                        hotspot_data[report][(report_nc, report_start, report_end)]['read_depth'][(int(position) - int(report_start))] = int(columns[_column_converter_snv['depth']])
                                    except KeyError:
                                        hotspot_data[report][(report_nc, report_start, report_end)]['read_depth'] = ["-"] * (int(report_end) - int(report_start) + 1)
                                        hotspot_data[report][(report_nc, report_start, report_end)]['read_depth'][(int(position) - int(report_start))] = int(columns[_column_converter_snv['depth']])
        try:
            for key, info in hotspot_data['hotspot'].items():
                variants_pindel = extract_hotspot_pindel_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped)
                if variants_pindel is not None:
                    for variant in variants_pindel:
                        output.write("\n" + "\t".join(variant))

                variants_bwa = extract_hotspot_snv_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped, '1-hotspot')
                if variants_bwa is not None:
                    for variant in variants_bwa:
                        output.write("\n" + "\t".join(variant))

                if variants_pindel is None and variants_bwa is None:
                    variants = no_variant_found(sample,info,key,read_depth_classes,annovar_header_mapper, ampliconmapped, '1-hotspot')
                    if variants is not None:
                        for variant in variants:
                            output.write("\n" + "\t".join(variant))
        except KeyError as e:
            if not e.args[0] == 'hotspot':
                raise e
            print("No hotspot data found")
        try:
            for key, info in hotspot_data['region'].items():
                variants_pindel = extract_hotspot_pindel_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped,'3-check')
                if variants_pindel is not None:
                    for variant in variants_pindel:
                        output.write("\n" + "\t".join(variant))

                variants_bwa = extract_hotspot_snv_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped, '3-check')
                if variants_bwa is not None:
                    for variant in variants_bwa:
                        output.write("\n" + "\t".join(variant))

        except KeyError as e:
            if not e.args[0] == 'region':
                raise e
            print("No region data found")


        try:
            for key, info in hotspot_data['indel'].items():
                variants_pindel = extract_hotspot_pindel_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped,'3-check')
                if variants_pindel is not None:
                    for variant in variants_pindel:
                        output.write("\n" + "\t".join(variant))

                variants_bwa = extract_hotspot_snv_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped, '3-check')
                if variants_bwa is not None:
                    for variant in variants_bwa:
                        output.write("\n" + "\t".join(variant))

        except KeyError as e:
            if not e.args[0] == 'indel':
                raise e
            print("No indel data found")
        try:
            for key, info in hotspot_data['region_all'].items():
                variants_pindel = extract_hotspot_pindel_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped,'3-check')
                if variants_pindel is not None:
                    for variant in variants_pindel:
                        output.write("\n" + "\t".join(variant))

                variants_bwa = extract_hotspot_snv_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped, '3-check')
                if variants_bwa is not None:
                    for variant in variants_bwa:
                        output.write("\n" + "\t".join(variant))

                #if variants_pindel is None and variants_bwa is None:
                variants = no_variant_found(sample,info,key,read_depth_classes,annovar_header_mapper, ampliconmapped,'3-check')
                if variants is not None:
                    for variant in variants:
                        output.write("\n" + "\t".join(variant))
        except KeyError as e:
            if not e.args[0] == 'region_all':
                raise e

def filtered_mutations_to_vcf(filtered_mutations_file, output_file, ampregion_file, nc2chr, min_vaf=0,print_all=False):
    vcf_header = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "##reference=hg19",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        ]
        )
    nc2chr = create_chr_mapper(nc2chr, False)
    regions = import_ampregion_seq(ampregion_file)
    with open(filtered_mutations_file, 'r') as input_data:
        with open(output_file, 'w') as output_data:
            output_data.write(vcf_header)
            header_map = {c[1]: c[0] for c in enumerate(input_data.readline().lstrip("#").split("\t"))}
            for line in input_data:
                line = line.rstrip("\n").rstrip("\r\n")
                if not line.startswith("#") and not len(line) == 0 and not line.lower().startswith("Sample"):
                    columns = line.split("\t")
                    if print_all or (columns[header_map['Variant_allele_ratio']] != '-' and float(columns[header_map['Variant_allele_ratio']]) > min_vaf):
                        id = "NA"
                        qual = "."
                        info = "."
                        filter_settings = "."

                        dbSnp = columns[header_map['dbSNP_id']]
                        cosmic = columns[header_map['Cosmic']]
                        if dbSnp is "-":
                            if cosmic is "-":
                                id = "."
                            else:
                                cosmic = cosmic.split(";")
                                cosmic = cosmic[0].split("=")
                                id = cosmic[1]
                        else:
                            id = dbSnp
                            if not cosmic is "-":
                                cosmic = cosmic.split(";")
                                cosmic = cosmic[0].split("=")
                                id = id + ";" + cosmic[1]
                        columns[header_map['Start']] = int(columns[header_map['Start']])
                        columns[header_map['End']] = int(columns[header_map['End']])
                        ref_seq = ""
                        var_seq = ""
                        if columns[header_map['Reference_base']] == "-":
                            try:
                                for s in regions[columns[header_map['Chr']]].keys():
                                    if s <= columns[header_map['Start']]:
                                        for e in regions[columns[header_map['Chr']]][s].keys():
                                            if e >= columns[header_map['End']]:
                                                rSeq = regions[columns[header_map['Chr']]][s][e][(columns[header_map['Start']] - s):(columns[header_map['Start']] - s + 1)]
                                                vSeq = rSeq + columns[header_map['Variant_base']]

                            except KeyError:
                                pass
                        elif columns[header_map['Variant_base']] == "-":
                            columns[header_map['Start']] -= 1
                            try:
                                for s in regions[columns[header_map['Chr']]].keys():
                                    if s <= columns[header_map['Start']]:
                                        for e in regions[columns[header_map['Chr']]][s].keys():
                                            if e >= columns[header_map['End']]:
                                                vSeq = rSeq = regions[columns[header_map['Chr']]][s][e][(columns[header_map['Start']] - s):(columns[header_map['Start']] - s + 1)]
                                                rSeq = rSeq + columns[header_map['Reference_base']]

                            except KeyError:
                                pass
                        else:
                            rSeq = columns[header_map['Reference_base']]
                            vSeq = columns[header_map['Variant_base']]

                        output_data.write("\n" + "\t".join([
                            nc2chr[columns[header_map['Chr']]],
                            str(columns[header_map['Start']]),
                            id,
                            rSeq,
                            vSeq,
                            qual,
                            filter_settings,
                            info,
                            "GT:DP:AD\t0/1:" + columns[header_map['Total_read_depth']] + ":" + columns[header_map['Reference_read_depth']] + "," + columns[header_map['Variant_read_depth']]]))


def calc_amplication(sample, tumourType, output, snpmania_variations_file, amplification_regions_file, background_regions_file, chr_to_nc):
    def add_data(data,chromosome,pos,rd):
        if chromosome in data:
            for start in data[chromosome]:
                if pos > start:
                    for end in data[chromosome][start]:
                        if pos <= end:
                            data[chromosome][start][end]['totDepth'] += rd
                            data[chromosome][start][end]['covBases'] += 1
        return data
    chr_to_nc_data = create_chr_mapper(chr_to_nc)
    amplification_regions = import_bed_file(amplification_regions_file, chr_to_nc_data, region_data_structure_generator)
    background_regions = import_bed_file(background_regions_file, chr_to_nc_data, region_data_structure_generator)
    with open(snpmania_variations_file, 'r') as snpmania_var_data:
        for line in snpmania_var_data:
            line = line.rstrip()
            if not line.startswith("#") and not line.startswith("\n"):
                # Populate referebce dictionary
                columns = line.split("\t")
                pos = columns[_column_converter_snv['position']]#int(columns[_column_converter_snv['position']])
                chromosome = columns[_column_converter_snv['nc_number']]
                depth = int(columns[_column_converter_snv['depth']])
                if chromosome in amplification_regions:
                    amplification_regions = add_data(amplification_regions, chromosome, pos, depth)
                else:
                    background_regions = add_data(background_regions, chromosome, pos, depth)
    import pprint
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(background_regions)
    with open(output,'w') as output:
        backTotDepth = 0
        backCovBases = 0

        for chromosome in background_regions:
            for start in background_regions[chromosome]:
                for end, data in background_regions[chromosome][start].items():
                    backTotDepth += data['totDepth']
                    backCovBases += data['covBases']
        backRdBase = 0
        if float(backCovBases > 0):
            backRdBase = float(backTotDepth) / float(backCovBases)
        for chromosome in amplification_regions:
            for start in amplification_regions[chromosome]:
                for end, data in amplification_regions[chromosome][start].items():
                    ampRdBase = 0
                    if data['covBases'] > 0:
                        ampRdBase = float(data['totDepth']) / float (data['covBases'])
                    amplification = 0
                    if backRdBase > 0:
                        amplification = float(ampRdBase) / float(backRdBase)
                    else:
                        amplification = "NA"
                    output.write("\t".join([
                        sample,
                        tumourType,
                        chromosome,
                        str(start),
                        str(end),
                        data['gene'],
                        str(amplification)]) + "\n")


if __name__ == "__main__":
    import doctest
    import logging
    import sys
    logging.basicConfig(level=logging.DEBUG, stream=sys.stdout, format='%(message)s')
    doctest.testmod()
    #result = doctest.testmod()
    #sys.exit(result.failed)
