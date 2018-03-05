# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import unittest
import shutil
import tempfile

class TestHandlersMethods(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()
        with open(os.path.join(self.tempdir,"mapper.tsv"),'w') as mapper:
            mapper.write("#Chr name\tNC\tID\tLength\n")
            mapper.write("chr1\tNC_000001.9\tChr1#NC_000001.9#1#247249719#-1\t247249719\n")
            mapper.write("chr2\tNC_000002.10\tChr2#NC_000002.10#1#242951149#-1\t242951149")

    def tearDown(self):
        # delete fixtures
        shutil.rmtree(self.tempdir)
        pass

    def test_chr_mapping(self):
        from scripts.lib.common.files.handlers import chr_mapping
        mapping = chr_mapping(os.path.join(self.tempdir, 'mapper.tsv'))
        print(mapping)
        self.assertTrue(len(mapping) == 2)
        self.assertTrue('NC_000001.9' in mapping.keys() and mapping['NC_000001.9'] == '1')
        self.assertTrue('NC_000002.10' in mapping.keys() and mapping['NC_000002.10'] == '2')

if __name__ == '__main__':
    unittest.main()
