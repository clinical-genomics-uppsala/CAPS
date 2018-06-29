#ToDo remove version found in wp1 parser
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

_column_bed_file = {'chr': 0, 'start': 1, 'end': 2, 'name': 3, 'score': 4, 'strand': 5}

def region_data_structure_generator(data, columns, mapper,start_position_modifier=1):
    chrom = columns[mapper['chr']]
    start = str(int(columns[mapper['start']]) + start_position_modifier) #int(columns[mapper['start']]) + start_position_modifier
    end = columns[mapper['end']] #int(columns[mapper['end']])
    name = columns[mapper['name']]
    if not chrom in data:
        data[chrom] = {}
        data[chrom][start] = {}
        data[chrom][start][end] = {}
        data[chrom][start][end]['gene'] = name
        data[chrom][start][end]['totDepth'] = 0
        data[chrom][start][end]['covBases'] = 0
    else:
        if not start in data[chrom]:
            data[chrom][start] = {}
            data[chrom][start][end] = {}
            data[chrom][start][end]['gene'] = name
            data[chrom][start][end]['totDepth'] = 0
            data[chrom][start][end]['covBases'] = 0
        else:
            if not end in data[chrom][start]:
                data[chrom][start][end] = {}
                data[chrom][start][end]['gene'] = name
                data[chrom][start][end]['totDepth'] = 0
                data[chrom][start][end]['covBases'] = 0
    return data

def import_bed_file(input_file, chr_to_nc, region_data_structure_generator):
    data = dict()
    with open(input_file, 'r') as lines:
        for line in lines:
            if line.startswith("chr"):
                columns = line.strip('\n').rstrip('\r').split("\t")
                columns[_column_bed_file['chr']] = chr_to_nc[columns[_column_bed_file['chr']].replace("chr","")]
                data = region_data_structure_generator(data,columns,_column_bed_file)
    return data


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
