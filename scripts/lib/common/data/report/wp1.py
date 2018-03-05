from scripts.lib.common.data.parser.wp1 import read_hotspot_file, read_blacklist, read_maintranscripts, read_multibpvariants, parse_reference_info, get_transcript_info,get_amplicon_info, extract_hotspot_pindel_info, extract_hotspot_snv_info
from scripts.lib.common.data.parser.wp1 import valid_amplicon_nonhotspot, valid_amplicon_hotspot, get_read_level, merge_comment, no_variant_found
from scripts.lib.common.data.parser.jsnpmania import _column_converter_snv
from scripts.lib.common.data.parser.annovar import generate_headermap
from scripts.lib.common.data.filters.annovar import depth_and_vaf, overlap_region, is_exonic
from scripts.lib.common.data.filters.annovar import is_splicing, is_nonsynonymous, in_1000g
from scripts.lib.common.data.filters.annovar import check_key_exists,contains_valid_amplicons,check_num_valid_ammplicons
from scripts.lib.common.data.filters.annovar import is_pindel_line
from scripts.lib.common.data.filters.wp1 import and_condition, or_condition, evalute_else, gene_name, insertion_var, vaf_above_or_equal, base_in_ref_or_var, evaluate_tuple

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
        ## print("Filter 1")
        if overlap_region(columns,intronic_data,headermap) or \
            (is_nonsynonymous(columns, headermap) and \
            (is_exonic(columns, headermap) or is_splicing(columns, headermap))):
            # Exclude, if it's found in the blacklist
            ## print("Filter 2")
            #if columns[headermap['Start']] == "29451783":
            #    print(str(key) + " " + str(key in blacklist_data))
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
    #TGFBR2_filter = (and_condition, (gene_name("TGFBR2"), (evalute_else, (insertion_var(1), vaf_above_or_equal(tgfbr2_1bp), vaf_above_or_equal(tgfbr2_2bp)))))
    #ACVR2A_filter = (and_condition, (gene_name("ACVR2A"), (evalute_else, (insertion_var(1), vaf_above_or_equal(acvr2a_1bp), vaf_above_or_equal(acvr2a_2bp)))))
    #final_filter = (and_condition, base_in_ref_or_var("A"), (or_condition, TGFBR2_filter, ACVR2A_filter))

    with open(output_file,'w') as output:
        output.write("#" + create_filtered_mutations_header())
        with open(input_file, 'r') as mutations:
            mapper = generate_headermap(mutations.readline(),"#Sample")
            for line in mutations:
                columns = line.rstrip("\n").rstrip("\r").split("\t")
                if filter_rows(columns, mapper):
                    output.write("\n" + line)


def generate_filtered_mutations(sample, output_file, hotspot_file, snpmania_variants_file, annovar_file, multibp_file, chr_to_nc_file, min_read_depth, min_vaf, read_depth_classes, genome1000_vaf, ampliconmapped, blacklist_file, transcript_file):
    """
        >>> hotspot = '/snakemake-workflows/resources/Mutations_Lung_20171219.csv'
        >>> snpmania = '/snakemake-workflows/resources/18-215.ampliconmapped.variations'
        >>> pindel = '/snakemake-workflows/resources/annovarOutput'
        >>> multibp = '/snakemake-workflows/resources/MultipleBp_variations.csv'
        >>> min_read_depth = 30
        >>> min_vaf = 0.01
        >>> read_depth_classes = [(300 ,"ok","yes"), (30 ,"low","yes"),(0 ,"low","not analyzable")]
        >>> max_1000genome = 0.02
        >>> blacklist = '/snakemake-workflows/resources/blacklist_variantClusterAllvarPindel_20151110.txt'
        >>> transcript = '/snakemake-workflows/resources/mainTranscripts.txt'
        >>> generate_filtered_mutations(hotspot, snpmania, pindel, multibp, min_read_depth, min_vaf, read_depth_classes, max_1000genome, True, blacklist, transcript)
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

                    if is_pindel_line(columns, annovar_header_mapper):
                        data_source = 'pindel'
                    for (amplicon_validator, reports) in [(valid_amplicon_hotspot, ('hotspot',)), (valid_amplicon_nonhotspot ,('region',))]:
                        if not added:
                            if data_source == 'pindel' or contains_valid_amplicons(amplicon_validator,columns, annovar_header_mapper):
                                for report in reports:
                                    if not added:
                                        if key in hotspot_data[report]:
                                            try:
                                                hotspot_data[report][key].get(data_source).append(columns)
                                            except KeyError as e:
                                                hotspot_data[report][key][data_source] = [columns]
                                            added = True
                                            break
                                        else:
                                            (chrom, start, end) = (columns[annovar_header_mapper['Chr']], columns[annovar_header_mapper['Start']], columns[annovar_header_mapper['End']])
                                            for entry_key in hotspot_data[report]:
                                                (entry_chrom, entry_start, entry_end) = entry_key
                                                if entry_chrom == chrom and \
                                                    start <= entry_end and end > entry_start:
                                                    try:
                                                        hotspot_data[report][entry_key].get(data_source).append(columns)
                                                    except:
                                                        hotspot_data[report][entry_key][data_source] = [columns]
                                                    added = True
                                                    break
                    if not added:
                        (ref_plus, ref_minus, var_plus, var_minus, ref_all, var_all) = get_amplicon_info(columns, annovar_header_mapper, ampliconmapped)

                        if data_source == 'pindel' or (contains_valid_amplicons(valid_amplicon_nonhotspot, columns, annovar_header_mapper) and ampliconmapped):
                            (found, level) = get_read_level(read_depth_classes, columns[annovar_header_mapper['Read_depth']])
                            (aa, cds, acc_num, exon, exonic_type, comment) =  get_transcript_info(columns, annovar_header_mapper, transcript_data)
                            comment = merge_comment(comment,None)
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

                variants_bwa = extract_hotspot_snv_info(info, read_depth_classes, annovar_header_mapper, transcript_data, ampliconmapped)
                if variants_bwa is not None:
                    for variant in variants_bwa:
                        output.write("\n" + "\t".join(variant))

                if variants_pindel is None and variants_bwa is None:
                    variants = no_variant_found(sample,info,key,read_depth_classes,annovar_header_mapper, '1-hotspot')
                    for variant in variants:
                        output.write("\n" + "\t".join(variant))
        except KeyError as e:
            if str(e) != 'hotspot':
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
            if str(e) != 'region':
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
            if str(e) != 'indel':
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

                if variants_pindel is None and variants_bwa is None:
                    variants = no_variant_found(sample,info,key,read_depth_classes,annovar_header_mapper,'3-check')
                    for variant in variants:
                        output.write("\n" + "\t".join(variant))
        except KeyError as e:
            if str(e) != 'region_all':
                raise e





if __name__ == "__main__":
    import doctest
    import logging
    import sys
    logging.basicConfig(level=logging.DEBUG, stream=sys.stdout, format='%(message)s')
    doctest.testmod()
    #result = doctest.testmod()
    #sys.exit(result.failed)
