# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import unittest
import shutil
import tempfile
import logging

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

class TestPindelMethods(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()
        with open(os.path.join(self.tempdir,"vcf"),'w') as variations:
            variations.write("##fileformat=VCFv4.0\n")
            variations.write("##fileDate=20090401\n")
            variations.write("##source=pindel\n")
            variations.write("##reference=Hg19\n")
            variations.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
            variations.write("##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">\n")
            variations.write("##INFO=<ID=PF,Number=1,Type=Integer,Description=\"The number of samples carry the variant\">\n")
            variations.write("##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">\n")
            variations.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
            variations.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
            variations.write("##INFO=<ID=NTLEN,Number=.,Type=Integer,Description=\"Number of bases inserted in place of deleted code\">\n")
            variations.write("##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
            variations.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            variations.write("##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference depth, how many reads support the reference\">\n")
            variations.write("##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allele depth, how many reads support this allele\">\n")
            variations.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
            variations.write("chr2\t140924351\t.\tA\tAAAAC\t.\tPASS\tEND=140924351;HOMLEN=23;HOMSEQ=AAACAAACAAACAAACAAACAAA;SVLEN=4;SVTYPE=INS\tGT:AD\t1/1:0,32\n")
            variations.write("chr3\t47058950\t.\tACTGCTGTCATGTAAGGTACGCATCCCTCCCCAAACCTTCCCTCCCCGTTC\tA\t.\tPASS\tEND=47059000;HOMLEN=3;HOMSEQ=CTG;SVLEN=-50;SVTYPE=DEL\tGT:AD\t0/0:1153,3\n")
            variations.write("chr19\t11095086\t.\tCTGCCCACTAGGGCTGCAGGCAGCCTCTGGACCGAGGGCCTTACTTGGAGGATGGGGGGAAGCCTTCTTGTTGGAGGTGTCCTGCCTTGGCT\tCA\t.\tPASS\tEND=11095177;HOMLEN=0;SVLEN=-91;SVTYPE=RPL;NTLEN=1\tGT:AD\t0/0:10805,2284\n")
        with open(os.path.join(self.tempdir,"deletions"),'w') as deletions:
            deletions.write("3455\tD 91\tNT 1 \"A\"\tChrID chr19\tBP 11095086\t11095178\tBP_range 11095086\t11095178\tSupports 2284\t11\t+ 1134\t5\t- 1150\t6\tS1 1306385\tSUM_MS 12514\t1\tNumSupSamples 1\t1\tsample1 10805 6496 1134 5 1150 6\n")
            deletions.write("737\tD 50\tNT 0 \"\"\tChrID chr3\tBP 47058950\t47059001\tBP_range 47058950\t47059004\tSupports 3\t1\t+ 3\t1\t- 0\t0\tS1 4\tSUM_MS 180\t1\tNumSupSamples 1\t1\tsample1 45 1153 3 1 0 0")
        with open(os.path.join(self.tempdir,"insertions"),'w') as insertions:
            insertions.write("105\tI 4\tNT 4 \"AAAC\"\tChrID chr2\tBP 140924351\t140924352\tBP_range 140924351\t140924375\tSupports 32\t2\t+ 15\t1\t- 17\t1\tS1 288\tSUM_MS 1920\t1\tNumSupSamples 1\t1\tsample1 0 0 15 1 17 1")

    def tearDown(self):
        # delete fixtures
        shutil.rmtree(self.tempdir)
        pass

    def test_ref_import(self):


        from scripts.lib.common.data.parser.pindel import convert_to_annovar_input
        convert_to_annovar_input(
            "sample1",
            os.path.join(self.tempdir, 'output'),
            os.path.join(self.tempdir, 'vcf'),
            os.path.join(self.tempdir, 'deletions'),
            os.path.join(self.tempdir, 'insertions'),
            20,
            0.01)
        self.maxDiff = None
        with open(os.path.join(self.tempdir, 'output'),'r') as result:
            line = result.readline()
            self.assertEqual(line,"2\t140924351\t140924351\t-\tAAAC\tcomments: sample=sample1 variantAlleleRatio=1.0 alleleFreq=0,32 readDepth=32 Tumor_Ins=+15|-17 Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=- Tumor_ref_minusAmplicons=-\n")
            line  = result.readline()
            self.assertEqual(line,"3\t47058951\t47059000\tCTGCTGTCATGTAAGGTACGCATCCCTCCCCAAACCTTCCCTCCCCGTTC\t-\tcomments: sample=sample1 variantAlleleRatio=0.0625 alleleFreq=45,3 readDepth=48 Tumor_Del=+3|-0 Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=- Tumor_ref_minusAmplicons=-\n")
            line = result.readline()
            self.assertEqual(line,"19\t11095087\t11095177\tTGCCCACTAGGGCTGCAGGCAGCCTCTGGACCGAGGGCCTTACTTGGAGGATGGGGGGAAGCCTTCTTGTTGGAGGTGTCCTGCCTTGGCT\tA\tcomments: sample=sample1 variantAlleleRatio=0.26013667425968107 alleleFreq=6496,2284 readDepth=8780 Tumor_Del=+1134|-1150 Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=- Tumor_ref_minusAmplicons=-\n")

if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    unittest.main()
