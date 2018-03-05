# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import unittest
import logging

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

class TestWP1ReportsMethods(unittest.TestCase):
    def setUp(self):
        self.filter1 = lambda column, mapper: True
        self.filter2 = filter2 = lambda column, mapper: False
        self.mapper = {'Gene': 0,'Reference_base': 1,'Variant_base':2,'Variant_allele_ratio':3,'Start':4,'Stop':5}
        self.columns = ["BRCA1","A","G",0.1,5,5]
        self.columns2 = ["BRCA1","A","GA",0.1,6,6]

    def tearDown(self):
        pass

    def test_evaluate_tuple(self):
        from scripts.lib.common.data.filters.wp1 import evaluate_tuple
        self.assertEqual(evaluate_tuple([],[],self.filter1), True)
        self.assertEqual(evaluate_tuple([],[],self.filter2), False)

    def test_evalute_else(self):
        from scripts.lib.common.data.filters.wp1 import evalute_else
        self.assertEqual(evalute_else(self.columns,self.mapper,(self.filter1, lambda c,m: True,lambda c,m: False)),True)
        self.assertEqual(evalute_else(self.columns,self.mapper,(self.filter2, lambda c,m: True,lambda c,m: False)),False)
        with self.assertRaises(Exception):
            evalute_else(self.columns,self.mapper,(filter1, lambda c,m: True))

    def test_and_condition(self):
        from scripts.lib.common.data.filters.wp1 import and_condition
        self.assertEqual(and_condition([],[],self.filter1,self.filter1), True)
        self.assertEqual(and_condition([],[],self.filter1,self.filter2), False)
        self.assertEqual(and_condition([],[],(and_condition,self.filter1,self.filter1),self.filter1), True)
        self.assertEqual(and_condition([],[],(and_condition,self.filter1,self.filter1),self.filter2), False)

    def test_or_condition(self):
        from scripts.lib.common.data.filters.wp1 import or_condition, and_condition
        self.assertEqual(or_condition([],[],self.filter1,self.filter1), True)
        self.assertEqual(or_condition([],[],self.filter1,self.filter2), True)
        self.assertEqual(or_condition([],[],self.filter2,self.filter2), False)
        self.assertEqual(or_condition([],[],(or_condition,self.filter2,self.filter2),self.filter1), True)
        self.assertEqual(or_condition([],[],(and_condition,self.filter1,self.filter2),self.filter2), False)

    def test_gene_name(self):
        from scripts.lib.common.data.filters.wp1 import gene_name
        self.assertEqual(gene_name("BRCA1")(self.columns,self.mapper), True)
        self.assertEqual(gene_name("BRCA2")(self.columns,self.mapper), False)

    def test_base_in_ref_or_var(self):
        from scripts.lib.common.data.filters.wp1 import base_in_ref_or_var
        self.assertEqual(base_in_ref_or_var("A")(self.columns,self.mapper), True)
        self.assertEqual(base_in_ref_or_var("G")(self.columns,self.mapper), True)
        self.assertEqual(base_in_ref_or_var("T")(self.columns,self.mapper), False)

    def test_vaf_above_or_equal(self):
        from scripts.lib.common.data.filters.wp1 import vaf_above_or_equal
        self.assertEqual(vaf_above_or_equal(0.1)(self.columns,self.mapper), True)
        self.assertEqual(vaf_above_or_equal(0.11)(self.columns,self.mapper), False)

    def test_insertion_var(self):
        from scripts.lib.common.data.filters.wp1 import insertion_var
        self.assertEqual(insertion_var(1)(self.columns,self.mapper), True)
        self.assertEqual(insertion_var(2)(self.columns,self.mapper), False)
        self.assertEqual(insertion_var(2)(self.columns2,self.mapper), True)

if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.INFO, stream=sys.stdout, format='%(message)s')
    unittest.main()
