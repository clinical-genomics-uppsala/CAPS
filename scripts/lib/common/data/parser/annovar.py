def generate_header(tumor_vs_normal = False):
    header = [
        "Sample",
        "Chr",
        "Start",
        "End",
        "Reference_base",
        "Variant_base",
        "Gene",
        "Type",
        "Exonic_type",
        "Variant_allele_ratio",
        "#reference_alleles",
        "#_variant_alleles",
        "Read_depth",
        "Ratio_in_1000Genome",
        "dbSNP_id",
        "Clinically_flagged_dbSNP",
        "ESP_6500",
        "Cosmic",
        "ClinVar_CLNDBN",
        "ClinVar_CLINSIG"]
    if not tumor_vs_normal:
        header += [
        "Strands_A",
        "Strands_G",
        "Strands_C",
        "Strands_T",
        "Strands_Ins",
        "Strands_Del",
        "#variant_+_amplicons",
        "#variant_-_amplicons",
        "#reference_+_amplicons",
        "#reference_-_amplicons",
        "Variant_ampliconinfo"]
    else:
        header += [
        "Tumor_Strands_A",
        "Tumor_Strands_G",
        "Tumor_Strands_C",
        "Tumor_Strands_T",
        "Tumor_Strands_Ins",
        "Tumor_Strands_Del",
        "Normal_Strands_A",
        "Normal_Strands_G",
        "Normal_Strands_C",
        "Normal_Strands_T",
        "#Tumor_variant_+_amplicons",
        "#Tumor_variant_-_amplicons",
        "#Tumor_reference_+_amplicons",
        "#Tumor_reference_-_amplicons",
        "#Normal_reference_+_amplicons",
        "#Normal_reference_-_amplicons",
        "Tumor_Variant_ampliconinfo",
        "Tumor_Reference_ampliconinfo"]
    header += [
        "Reference_ampliconinfo",
        "Transcripts"]
    return header
"""
    >>> extract_gene_info("U2AF1","splicing","NM_006758:exon2:c.44+4A>G","-")
    ("U2AF1","U2AF1:NM_006758:exon2:c.44+4A>G")
"""

def extract_gene_info(gene, func_type, transcript, aa_change):
    if "splicing" in func_type:
        if "exonic;splicing" in func_type:
            if "ncRNA" in func_type:
                gene = gene.split(";")[0]
                transcript += "," + aa_change
            else:
                gene = gene.split(";")[1]
        return (gene, transcript.replace("NM",gene + ":NM"))
    else:
        return (gene, aa_change)

#ToDo move to more generic file
def generate_headermap(line,startswith="Chr", sep="\t"):
    """
        >>> line = "Chr\\tStart\\tEnd\\tRef\\tAlt\\tFunc.refGene\\tGene.refGene\\tGeneDetail.refGene\\tExonicFunc.refGene\\tAAChange.refGene\\tsnp138\\tsnp138NonFlagged\\tesp6500siv2_ea\\tcosmic70\\tclinvar_20150629\\tOtherinfo"
        >>> generate_headermap(line)
        {'Chr': 0, 'Start': 1, 'End': 2, 'Ref': 3, 'Alt': 4, 'Func.refGene': 5, 'Gene.refGene': 6, 'GeneDetail.refGene': 7, 'ExonicFunc.refGene': 8, 'AAChange.refGene': 9, 'snp138': 10, 'snp138NonFlagged': 11, 'esp6500siv2_ea': 12, 'cosmic70': 13, 'clinvar_20150629': 14, 'Otherinfo': 15}
    """
    if not line.startswith(startswith):
        raise Exception("Header line should start with \"{0}\"".format(startswith))
    else:
        if line.startswith("#"):
            line = line[1:]
        return dict([(v, i) for i,v in enumerate(line.rstrip().split(sep))])

def process_clinvar_field(data):
    """
        >>> data = ""
        >>> process_clinvar_field(data)
        {}
        >>> data = "CLINSIG=untested;CLNDBN=Breast-ovarian_cancer\x2c_familial_1;CLNREVSTAT=no_assertion_provided;CLNACC=RCV000112677.1;CLNDSDB=GeneReviews:MedGen:OMIM:Orphanet;CLNDSDBID=NBK1247:C2676676:604370:ORPHA145"
        >>> process_clinvar_field(data)
        {'CLINSIG': 'untested', 'CLNDBN': 'Breast-ovarian_cancer\x2c_familial_1'}
    """
    clinvar_data = {}
    if len(data) > 0:
        for d in data.split(";"):
            if "CLNDBN" in d or "CLINSIG" in d:
                clinvar_data.update(dict([tuple(d.split("="))]))
    return clinvar_data

def parse_multianno_line(line, headermap):
    """
        >>> headermap = {"Func.refGene": 5, "Gene.refGene": 6, "GeneDetail.refGene": 7,'AAChange.refGene': 9}
        >>> line = "21\\t44527557\\t44527557\\tT\\tC\\tsplicing\\tU2AF1\\tNM_006758:exon2:c.44+4A>G,NM_001025203:exon2:c.44+4A>G\\t-\\t-\\trs372157069\\trs372157069\\t-\\t-\\t-\\tcomments: sample=sample1 variantAlleleRatio=0.0117647058823529 alleleFreq=84,1 readDepth=85 Tumor_A=0|0|0|0 Tumor_G=0|0|0|0 Tumor_C=1|0|0|0 Tumor_T=71|0|13|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=chr21:44527551-44527782:+:1 Tumor_ref_ampliconInfo=chr21:44527547-44527778:-:13#chr21:44527551-44527782:+:65"
        >>> parse_multianno_line(line,headermap)
        ('U2AF1', 'U2AF1:NM_006758:exon2:c.44+4A>G,U2AF1:NM_001025203:exon2:c.44+4A>G', {'sample': 'sample1', 'variantAlleleRatio': '0.0117647058823529', 'alleleFreq': '84,1', 'readDepth': '85', 'Tumor_A': '0|0|0|0', 'Tumor_G': '0|0|0|0', 'Tumor_C': '1|0|0|0', 'Tumor_T': '71|0|13|0', 'Tumor_var_plusAmplicons': '0', 'Tumor_var_minusAmplicons': '0', 'Tumor_ref_plusAmplicons': '1', 'Tumor_ref_minusAmplicons': '1', 'Tumor_var_ampliconInfo': 'chr21:44527551-44527782:+:1', 'Tumor_ref_ampliconInfo': 'chr21:44527547-44527778:-:13#chr21:44527551-44527782:+:65'})
    """
    if not line.startswith("Chr"):
        line = line.rstrip().split("\t")
        # Split and remote comments part
        comments = line[-1].split(" ")[1:]
        info = dict(map(lambda x: tuple(x.split("=")),comments))
        (gene, transcript) = extract_gene_info(
                                line[headermap['Gene.refGene']],
                                line[headermap['Func.refGene']],
                                line[headermap['GeneDetail.refGene']],
                                line[headermap['AAChange.refGene']])
        line[headermap['Func.refGene']].replace("cRNA_exonic;","")
        return (line, gene, transcript, info)

#def parse_annovar_output_line(line, headermap):
#    """
#        >>> heade = "#Sample\\tChr\\tStart\\tEnd\\tReference_base\\tVariant_base\\tGene\\tType\\tExonic_type\\tVariant_allele_ratio\\t#reference_alleles\\t#_variant_alleles\\tRead_depth\\tRatio_in_1000Genome\\tdbSNP_id\\tClinically_flagged_dbSNP\\tESP_6500\\tCosmic\\tClinVar_CLNDBN\\tClinVar_CLINSIG S#trands_A\\tStrands_G\\tStrands_C\\tStrands_T\\tStrands_Ins\\tStrands_Del\\t#variant_+_amplicons\\t#variant_-_amplicons\\t#reference_+_amplicons\\t#reference_-_amplicons\\tVariant_ampliconinfo\\tReference_ampliconinfo\\tTranscripts"
#        >>> line = "21\\t44527557\\t44527557\\tT\\tC\\tsplicing\\tU2AF1\\tNM_006758:exon2:c.44+4A>G,NM_001025203:exon2:c.44+4A>G\\t-\\t-\\trs372157069\\trs372157069\\t-\\t-\\t-\\tcomments: sample=sample1 variantAlleleRatio=0.0117647058823529 alleleFreq=84,1 readDepth=85 Tumor_A=0|0|0|0 Tumor_G=0|0|0|0 Tumor_C=1|0|0|0 Tumor_T=71|0|13|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=chr21:44527551-44527782:+:1 T#umor_ref_ampliconInfo=chr21:44527547-44527778:-:13#chr21:44527551-44527782:+:65"
#        >>> parse_multianno_line(line,headermap)
#        ('U2AF1', 'U2AF1:NM_006758:exon2:c.44+4A>G,U2AF1:NM_001025203:exon2:c.44+4A>G', {'sample': 'sample1', 'variantAlleleRatio': '0.0117647058823529', 'alleleFreq': '84,1', 'readDepth': '85', 'Tumor_A': '0|0|0|0', 'Tumor_G': '0|0|0|0', 'Tumor_C': '1|0|0|0', 'Tumor_T': '71|0|13|0', 'Tumor_var_plusAmplicons': '0', 'Tumor_var_minusAmplicons': '0', 'Tumor_ref_plusAmplicons': '1', 'Tumor_ref_minusAmplicons': '1', 'Tumor_var_ampliconInfo': 'chr21:44527551-44527782:+:1', 'Tumor_ref_ampliconInfo': 'chr21:44527547-44527778:-:13#chr21:44527551-44527782:+:65'})
#    """
#    if not line.startswith("Chr"):
#        line = line.rstrip().split("\t")
#        # Split and remote comments part
#        comments = line[-1].split(" ")[1:]
#        info = dict(map(lambda x: tuple(x.split("=")),comments))
#        (gene, transcript) = extract_gene_info(
#                                line[headermap['Gene.refGene']],
#                                line[headermap['Func.refGene']],
#                                line[headermap['GeneDetail.refGene']],
#                                line[headermap['AAChange.refGene']])
#        line[headermap['Func.refGene']].replace("cRNA_exonic;","")
#        return (gene, transcript, info)


def process_annovar_multianno_file(output_file, multianno_file, tumor_vs_normal = False, amplicon_mapped = False):
    with open(multianno_file, 'r') as multianno_lines:
        with open(output_file, 'w') as output:
            header_line = multianno_lines.readline()
            headermap = generate_headermap(header_line)
            output.write("#"+ "\t".join(generate_header(tumor_vs_normal)))
            for line in multianno_lines:
                (line, gene, transcript, info) = parse_multianno_line(line,headermap)
                clinvar_data = process_clinvar_field(line[headermap['clinvar_20150629']])
                allele_freq = info.get('alleleFreq','-,').split(",")
                data = [info.get('sample','-'),
                line[headermap['Chr']],
                line[headermap['Start']],
                line[headermap['End']],
                line[headermap['Ref']],
                line[headermap['Alt']],
                gene,
                line[headermap['Func.refGene']],
                line[headermap['ExonicFunc.refGene']],
                info.get('variantAlleleRatio','-'),
                allele_freq[0],
                allele_freq[1],
                info.get('readDepth','-,'),
                line[headermap['1000g2015aug_eur']],
                line[headermap['snp138']],
                "Yes" if line[headermap['snp138']].startswith("rs") and "-" in line[headermap['snp138NonFlagged']] else "No",
                line[headermap['esp6500siv2_ea']],
                line[headermap['cosmic70']],
                clinvar_data.get('CLNDBN','-'),
                clinvar_data.get('CLINSIG','-'),
                info.get('Tumor_A','-'),
                info.get('Tumor_G','-'),
                info.get('Tumor_C','-'),
                info.get('Tumor_T','-'),
                info.get('Tumor_Ins','-'),
                info.get('Tumor_Del','-')]
                if tumor_vs_normal:
                    data += [
                        info.get('Normal_A','-'),
                        info.get('Normal_G','-'),
                        info.get('Normal_C','-'),
                        info.get('Normal_T','-')
                    ]
                data += [
                    info.get('Tumor_var_plusAmplicons','-'),
                    info.get('Tumor_var_minusAmplicons','-'),
                    info.get('Tumor_ref_plusAmplicons','-'),
                    info.get('Tumor_ref_minusAmplicons','-')]
                if tumor_vs_normal:
                    data += [
                        info.get('Normal_ref_plusAmplicons','-'),
                        info.get('Normal_ref_minusAmplicons','-')]
                data += [
                    info.get('Tumor_var_ampliconInfo','-'),
                    info.get('Tumor_ref_ampliconInfo','-')]
                if tumor_vs_normal:
                    data += [info.get('Tumor_var_ampliconInfo','-')]
                data += [transcript]
                output.write("\n" + "\t".join(data))


if __name__ == "__main__":
    import doctest
    #import logging
    #import sys
    #logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    doctest.testmod()
    #result = doctest.testmod()
    #sys.exit(result.failed)
