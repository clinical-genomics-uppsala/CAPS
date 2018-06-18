# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import unittest
import shutil
import tempfile
import logging

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

class TestWP1ReportsMethods(unittest.TestCase):
    def setUp(self):
        from scripts.lib.common.data.report.wp1 import create_filtered_mutations_header
        self.tempdir = tempfile.mkdtemp()
        self.line1 = "18-99\tTGFBR2\tframeshift deletion\texon3\tp.E125fs\tc.374delA\tNM_003242\t-\t4-other\tyes\tok\t4328\t3039\t1273\t0.29413123844732\t-\t-\t-\tNo\tID=COSM1180952,COSM1744957;OCCURENCE=1(biliary_tract),1(large_intestine)\t-\t-\t4\t3\t3\t3\t-\t-\t-\t-\t-\t403|250|244|392\tchr3:30691813-30691914:-:684#chr3:30691813-30691914:+:578#chr3:30691860-30691926:+:445#chr3:30691860-30691926:-:501#chr3:30691492-30691945:+:13#chr3:30691496-30691948:-:22#chr3:30691691-30691968:+:4#chr3:30691795-30691909:+:773\t(0,0):chr3:30691813-30691914:-:288#(0,0):chr3:30691813-30691914:+:311#(0,0):chr3:30691860-30691926:+:140#(0,0):chr3:30691860-30691926:-:180#(0,0):chr3:30691492-30691945:+:4#(0,0):chr3:30691496-30691948:-:11#(0,0):chr3:30691795-30691909:+:332\tNC_000003.11\t30691872\t30691872\tA\t-\tTGFBR2:NM_003242:exon3:c.374delA:p.E125fs,TGFBR2:NM_001024847:exon4:c.449delA:p.E150fs"
        self.line2 = "18-99\tBRAF\tsplicing\t-\t-\t-\t-\t-\t4-other\tyes\tok\t2707\t1967\t51\t0.0188400443295161\t-\t-\t-\tNo\t-\t-\t-\t3\t4\t1\t2\t-\t-\t-\t-\t-\t99|3|637|1\tchr7:140434545-140434681:+:120#chr7:140434545-140434683:-:1659#chr7:140434529-140434775:-:49#chr7:140434525-140434771:+:5#chr7:140434552-140434838:+:3#chr7:140434552-140434838:-:11#chr7:140434552-140434895:-:10#chr7:140434427-140434694:+:25\t(0,1):chr7:140434545-140434681:+:7#(0,1):chr7:140434545-140434683:-:34#(0,1):chr7:140434529-140434775:-:6#(0,1):chr7:140434525-140434771:+:1#(0,1):chr7:140434552-140434838:-:1#(0,1):chr7:140434552-140434895:-:1\tNC_000007.13\t140434575\t140434576\tAA\t-\t-"

        with open(os.path.join(self.tempdir,"mutations"),'w') as mutations:
            mutations.write("#" + create_filtered_mutations_header())
            mutations.write("\n" + self.line1)
            mutations.write("\n" + self.line2)


    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_filter_mutations(self):
        from scripts.lib.common.data.report.wp1 import generate_filtered_mutations
        sample = "18-215"
        hotspot = '/snakemake/resources/Mutations_Lung_20171219.csv'
        snpmania = '/snakemake/resources/18-215.ampliconmapped.variations'
        pindel = ['/snakemake/resources/annovarOutput']
        multibp = '/snakemake/resources/MultipleBp_variations.csv'
        chr_to_nc = '/snakemake/resources/reference.chr_to_nc.hg19.info'
        output = '/snakemake/resources/output'
        min_read_depth = 30
        min_vaf = 0.01
        read_depth_classes = [(300 ,"ok","yes"), (30 ,"low","yes"),(0 ,"lowk","not analyzable")]
        max_1000genome = 0.02
        blacklist = '/snakemake/resources/blacklist_variantClusterAllvarPindel_20151110.txt'
        transcript = '/snakemake/resources/mainTranscripts.txt'
        generate_filtered_mutations(sample, output, hotspot, snpmania, pindel, multibp, chr_to_nc, min_read_depth, min_vaf, read_depth_classes, max_1000genome, True, blacklist, transcript)


    def test_extract_msi_markers(self):
        from scripts.lib.common.data.report.wp1 import extract_msi_markers
        extract_msi_markers(os.path.join(self.tempdir,"mutations"), os.path.join(self.tempdir,"msi"))

        with open(os.path.join(self.tempdir,"msi"),'r') as msi:
            from scripts.lib.common.data.report.wp1 import create_filtered_mutations_header
            self.maxDiff = None
            self.assertEqual(msi.readline().rstrip("\n").rstrip("\r"),"#"+create_filtered_mutations_header())
            self.assertEqual(msi.readline().rstrip("\n").rstrip("\r"),self.line1)
            self.assertEqual(msi.readline(),"")


if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.INFO, stream=sys.stdout, format='%(message)s')
    unittest.main()
