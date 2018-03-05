from scripts.lib.common.data.parser.jsnpmania import _extract_ref_variant_info, _extract_deletions, _extract_insertion
from scripts.lib.common.data.parser.pindel import parse_pindel_line
import vcf.filters
import collections

_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])

# ToDo remove sample_name, since it's only used to build up the dictionary
class ReadSupportAlternativeAlleleAnnotater:
    """
        Function used to add information about supporting reads for variants,
        supporting alleles forward (SAF) and supporting alleles reverse (SAR).
    """
    def __init__(self, sample_name, pindel_deletion_file , pindel_insertion_file):
        """
            :param sample_name: name of sample,
            :param pindel_deletion_file: deletion file from pindel, ".indels_D"
            :param pindel_insertion_file: insertion file from pindel, ".indels_SI"
        """
        self.deletion_data = dict() #used format { sample_name: {'chr1:1000:1002:-': {'fwd': 1, "rvd": 1}}}
        self.insertion_data = dict()#used format { sample_name: {'chr1:1000:1001:ACGT': {'fwd': 1, "rvd": 1}}}
        for sample in sample_name:
            self.deletion_data[sample] = dict()
            self.insertion_data[sample] = dict()
            with open(pindel_deletion_file[sample],'r') as deletions:
                for line in deletions:
                    if line and line[0].isdigit():
                        data = parse_pindel_line(line)
                        self.deletion_data[sample]["%s:%d:%d:%s" % (data['chr'],data['start']-1,data['stop'] - 1,data['indel'].replace("-",""))] = data
            with open(pindel_insertion_file[sample],'r') as insertions:
                for line in insertions:
                    if line and line[0].isdigit():
                        data = parse_pindel_line(line)
                        self.insertion_data[sample]["%s:%d:%d:%s" % (data['chr'],data['start']-1,data['stop']-1,data['indel'])] = data

    def annotate(self,input_file, output_file):
        """
            Create new vcf file containing tag SAF and SAR. Process the input
            vcf file and calculate SAF and SAR values for each record.

            :param input_file: vcf file that will be annotated
            :param input_file: new annotated vcf file
        """
        with open(input_file,'r') as input_vcf_file:
            with open(output_file,'w') as output_vcf_file:
                vcf_reader = vcf.Reader(input_vcf_file)
                vcf_reader.formats['SAF'] = _Format('SAF',1,'Integer',"Alternate allele observations on the forward strand")
                vcf_reader.formats['SAR'] = _Format('SAR',1,'Integer',"Alternate allele observations on the reverse strand")
                vcf_writer = vcf.Writer(output_vcf_file, vcf_reader)
                for record in vcf_reader:
                    self.annotate_record(record)
                    vcf_writer.write_record(record)

    def annotate_record(self,record):
        """
            Add the actual information to the record.

            :param record: the record that will be updated.
        """
        record.add_format('SAF')
        record.add_format('SAR')
        call_data = collections.namedtuple('calldata', record.FORMAT.split(":"))
        for sample_index in range(len(record.samples)):
            for allele in record.alleles[1:]:
                data = self.deletion_data
                if record.is_sv and record.var_subtype == "INS":
                    data = self.insertion_data
                elif not (record.is_deletion or (record.is_sv and record.var_subtype in ["DEL", "RPL"])):
                    raise Exception() # Raise an exception if we have a record type that we haven't considered.
                # Retrieved the stored pindel information.
                indel = data[record.samples[sample_index].sample]["%s:%d:%d:%s" % (record.CHROM,record.start,record.end,str(allele)[1:])]
                # Add the forward and reverse information
                prev_data = [ d for d in record.samples[sample_index].data] + [indel['fwd'], indel['rve']]
                record.samples[sample_index].data = call_data(*prev_data)

#from scripts.lib.common.data.annotaters.vcf import AlternativeAlleleReadSupportAnnotater; annotater = AlternativeAlleleReadSupportAnnotater('data/PindelOut/17-1920.indels_D', 'data/PindelOut/17-1920.indels_SI'); annotater.annotate('data/PindelAnnovar/17-1920.pindel.vcf','data/PindelAnnovar/17-1920.pindel.v2.vcf')
# def add_amplicon_information(vcf, sample_name, nc_to_chr, jsnpmania_variants, jsnpmania_insertions, jsnpmania_deletions,amplicon_min_depth):
#     reference = dict()
#     insertions = dict()
#     deletions = dict()
#     for keys in nc_to_chr:
#         reference[keys] = {}
#     with open(jsnpmania_variants,'r') as j_variants:
#         for line in j_variants:
#             line = line.rstrip()
#             if not line.startswith("#") and not line.startswith("\n") != 0:
#                 reference = _extract_ref_variant_info(reference, line, amplicon_min_depth)
#
#     with open(jsnpmania_insertions,'r') as j_insertions:
#         for line in j_insertions:
#             line = line.rstrip()
#             if not line.startswith("#") and not line.startswith("\n") != 0:
#                 alleles = _extract_insertion(sample_name, reference, line, nc_to_chr, amplicon_min_depth)
#                 for key, info in alleles.items():
#                     (nc, start, stop, ref, var) = key.split("#")
#                     columns = info.split("\t")[19].split()
#                     upper
#
#                     try:
#                         if int(start) in variants[nc]:
#                             insertions[nc_to_chr[chr]][int(start)][int(stop)] = var
#                         else:
#                             insertions[nc_to_chr[nc]][int(start)] = {int(stop): info}
#                     except KeyError:
#                         insertions[nc_to_chr[nc]] = {int(start): {int(stop): info}}
#     with open(jsnpmania_deletions,'r') as j_deletions:
#         for line in j_deletions:
#             line = line.rstrip()
#             if not line.startswith("#") and not line.startswith("\n") != 0:
#                 alleles = _extract_deletions(sample_name, reference, line, nc_to_chr, amplicon_min_depth)
#                 for key, info in alleles.items():
#                     (nc, start, stop, ref, var) = key.split("#")
#                     #try:
#                     #    variants[nc][pos][key] = info
#                     try:
#                         if int(start) in deletions[nc]:
#                             deletions[nc_to_chr[nc]][int(start)][int(stop)] = var
#                         else:
#                             deletions[nc_to_chr[nc]][int(start)] = {int(stop): info}
#                     except KeyError:
#                         deletions[nc_to_chr[nc]] = {int(start): {int(stop): info}}
