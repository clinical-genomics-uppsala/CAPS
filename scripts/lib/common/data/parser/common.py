def create_chr_mapper(mapper_file,  chr_to_nc=True):
    """
        First and second column of input file should be of the following
        format:
        First column: chr1
        Second column: NC_000001.10

        :param file: path to input file
        :return dict= {'NC_000001.10': 'chr1'} or {'chr1': 'NC_000001.10'}

    """
    chr_mapper = dict()
    with open(mapper_file,'r') as lines:
        for line in lines:
            if not line.startswith("#"):
                line = line.rstrip()
                columns = line.split("\t")
                if chr_to_nc:
                    chr_mapper[columns[0]] = columns[1]
                else:
                    chr_mapper[columns[1]] = columns[0]
    return chr_mapper

_column_converter_ampregion = {
    'chr': 1,
    'start': 2,
    'end': 3,
    'seq': 4}

def import_ampregion_seq(ampregion_file):
    """
        Convert ampregion file to dict.

        :param file: path to  tab separated input file with format name, chr,
                        start, end, seq
        return dict= {'NC_000001.10': {'1000': {'1050': 'ACGT...AGCT'}}}

    """
    regions = dict()
    with open(ampregion_file, 'r') as lines:
        for line in lines:
            columns = line.rstrip('\r\n').rstrip('\n').split("\t")
            if not columns[_column_converter_ampregion['chr']].startswith("NC_"):
                raise Exception("Chromosome column should use NC format, ex NC_000001.10")

            chrom = columns[_column_converter_ampregion['chr']]
            start = int(columns[_column_converter_ampregion['start']])
            end = int(columns[_column_converter_ampregion['end']])
            seq = columns[_column_converter_ampregion['seq']]

            if chrom in regions:
                if start in regions[chrom]:  # Check if start exists in the dictionary
                    if end > regions[chrom][start]:  # Check if the end is greater than the one already existing
                        regions[chrom][start] = {key: value for key, value in regions[chrom][start].items() if value != regions[chrom][start]}
                        regions[chrom][start][end] = seq
                else:  # If start doesn't exist in dict create it
                    regions[chrom][start] = {}
                    regions[chrom][start][end] = seq  # save sequence
            else:  # If chrom doesn't exist in dict create the whole structure
                regions[chrom] = {}
                regions[chrom][start] = {}
                regions[chrom][start][end] = seq
    return regions
