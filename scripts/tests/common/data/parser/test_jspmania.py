# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import unittest
import shutil
import tempfile
import logging

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

class TestJSNPManiaMethods(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()
        with open(os.path.join(self.tempdir,"variations"),'w') as variations:
            variations.write("6463\t3532\t0.5464954355562432\tNC_000014.8\t105258892\tC\t4|0|3532|2927\tC\tT\t0.5468338752128813\t0|2|2|0\t0|0|0|0\t1077|691|695|1069\t842|627|629|829\tchr14:105258855-105258982:-:2#chr14:105258810-105258944:-:2|0|chr14:105258859-105258986:+:1544#chr14:105258855-105258982:-:1010#chr14:105258810-105258944:-:324#chr14:105258806-105258944:+:363#chr14:105258778-105259019:-:52#chr14:105258774-105259016:+:228|chr14:105258859-105258986:+:1193#chr14:105258855-105258982:-:916#chr14:105258810-105258944:-:274#chr14:105258806-105258944:+:323#chr14:105258778-105259019:-:47#chr14:105258774-105259016:+:145\n")
            variations.write("6455\t1\t1.5491866769945777E-4\tNC_000014.8\t105258893\tA\t1|6454|0|0\tG\tA\t0.9998450813323005\t0|0|1|0\t1915|1318|1320|1901\t0|0|0|0\t0|0|0|0\tchr14:105258855-105258982:-:1|chr14:105258859-105258986:+:2737#chr14:105258855-105258982:-:1929#chr14:105258810-105258944:-:595#chr14:105258806-105258944:+:686#chr14:105258778-105259019:-:95#chr14:105258774-105259016:+:374|0|0\n")
            variations.write("6434\t6434\t1.0\tNC_000014.8\t105258894\tT\t0|0|0|6434\tT\tA/C/G\t1.0\t0|0|0|0\t0|0|0|0\t0|0|0|0\t1919|1311|1316|1888\t0|0|0|chr14:105258859-105258986:+:2731#chr14:105258855-105258982:-:1917#chr14:105258810-105258944:-:595#chr14:105258806-105258944:+:689#chr14:105258778-105259019:-:96#chr14:105258774-105259016:+:368\n")
            variations.write("78\t76\t0.9743589743589743\tNC_000001.10\t115252473\tT\t1|0|1|76\tT\tA/C\t0.987012987012987\t0|0|1|0\t0|0|0|0\t1|0|0|0\t43|0|33|0\tchr1:115252332-115252648:-:1|0|chr1:115252332-115252647:+:1|chr1:115252332-115252647:+:43#chr1:115252332-115252648:-:33\n")
            variations.write("1706\t1601\t0.9384525205158265\tNC_000010.10\t123256344\tT\t104|1|0|1601\tT\tA\t0.9390029325513196\t0|52|52|0\t0|0|0|1\t0|0|0|0\t235|567|568|231\t0|chr10:123256234-123256344:+:1|0|chr10:123256234-123256344:+:416#chr10:123256234-123256346:-:588\n")
            variations.write("1148\t795\t0.6925087108013938\tNC_000010.10\t123256345\tT\t353|0|0|795\tT\tA\t0.6925087108013938\t0|182|171|0\t0|0|0|0\t0|0|0|0\t23|374|374|24\t0|0|0|chr10:123256234-123256344:+:1#chr10:123256234-123256346:-:587\n")
            variations.write("670\t567\t0.8462686567164179\tNC_000010.10\t123256346\tA\t567|0|0|103\tA\tT\t0.8462686567164179\t1|292|274|0\t0|0|0|0\t0|0|0|0\t0|57|45|1\tchr10:123256234-123256346:-:564|0|0|chr10:123256234-123256346:-:3\n")
            variations.write("1412\t700\t0.5042492917847026\tNC_000010.10\t123239112\tG\t712|700|0|0\tG\tA\t0.5042492917847026\t170|83|19|428\t213|72|25|402\t0|0|0|0\t0|0|0|0\tchr10:123239059-123239249:-:36#chr10:123239055-123239245:+:340#chr10:123238910-123239122:+:258#chr10:123238910-123239121:-:64|chr10:123239059-123239249:-:40#chr10:123239055-123239245:+:406#chr10:123238910-123239122:+:192#chr10:123238910-123239121:-:47|0|0\n")
            variations.write("91\t89\t0.978021978021978\tNC_000007.13\t140499713\tA\t89|0|0|0\tA\tC/G/T\t1.0\t50|0|33|6\t0|0|0|0\t0|0|0|0\t0|0|0|0\tchr7:140499687-140500177:-:33#chr7:140499687-140500177:+:44|0|0|0\n")
            variations.write("227\t216\t0.9515418502202643\tNC_000007.13\t140508149\tT\t4|2|0|216\tT\tA\t0.9818181818181818\t0|4|0|0\t0|2|0|0\t0|0|0|0\t36|21|159|0\tchr7:140508026-140508255:-:4|chr7:140508026-140508255:-:2|0|chr7:140508000-140508412:+:36#chr7:140508001-140508412:-:19#chr7:140508026-140508255:-:155\n")
            variations.write("232\t225\t0.9698275862068966\tNC_000007.13\t140508150\tC\t0|0|225|0\tC\tA/G/T\t1.0\t0|0|0|0\t0|0|0|0\t34|43|148|0\t0|0|0|0\t0|0|chr7:140508000-140508412:+:34#chr7:140508001-140508412:-:18#chr7:140508026-140508255:-:166|0\n")
            variations.write("243\t221\t0.9094650205761317\tNC_000007.13\t140508151\tC\t7|0|221|0\tC\tA\t0.9692982456140351\t0|7|0|0\t0|0|0|0\t26|45|150|0\t0|0|0|0\tchr7:140508026-140508255:-:6|0|chr7:140508000-140508412:+:26#chr7:140508001-140508412:-:17#chr7:140508026-140508255:-:172|0\n")
            variations.write("265\t205\t0.7735849056603774\tNC_000007.13\t140508152\tA\t205|0|3|0\tA\tC\t0.9855769230769231\t0|45|160|0\t0|0|0|0\t0|3|0|0\t0|0|0|0\tchr7:140508001-140508412:-:19#chr7:140508026-140508255:-:177|0|chr7:140508026-140508255:-:3|0\n")
            variations.write("260\t217\t0.8346153846153846\tNC_000007.13\t140508153\tA\t217|0|1|0\tA\tC\t0.9954128440366973\t0|75|142|0\t0|0|0|0\t0|1|0|0\t0|0|0|0\tchr7:140508026-140508255:-:206|0|chr7:140508026-140508255:-:1|0\n")
            variations.write("262\t238\t0.9083969465648855\tNC_000007.13\t140508154\tA\t238|0|0|0\tA\tC/G/T\t1.0\t0|96|142|0\t0|0|0|0\t0|0|0|0\t0|0|0|0\tchr7:140508026-140508255:-:227|0|0|0\n")
            variations.write("273\t259\t0.9487179487179487\tNC_000007.13\t140508155\tA\t259|0|0|0\tA\tC/G/T\t1.0\t0|117|142|0\t0|0|0|0\t0|0|0|0\t0|0|0|0\tchr7:140508026-140508255:-:247|0|0|0\n")
            variations.write("275\t269\t0.9781818181818182\tNC_000007.13\t140508156\tA\t269|0|0|0\tA\tC/G/T\t1.0\t0|130|139|0\t0|0|0|0\t0|0|0|0\t0|0|0|0\tchr7:140508026-140508255:-:258|0|0|0\n")
            variations.write("281\t279\t0.9928825622775801\tNC_000007.13\t140508157\tA\t279|0|0|0\tA\tC/G/T\t1.0\t0|138|141|0\t0|0|0|0\t0|0|0|0\t0|0|0|0\tchr7:140508026-140508255:-:267|0|0|0\n")
            variations.write("4199\t1711\t0.40747797094546323\tNC_000004.11\t55968053\tA\t1711|0|2487|0\tC\tA\t0.5924249642686994\t428|426|426|431\t0|0|0|0\t654|590|575|668\t0|0|0|0\tchr4:55968023-55968093:+:384#chr4:55967974-55968185:+:243#chr4:55967975-55968186:-:188#chr4:55967939-55968092:+:232#chr4:55967938-55968091:-:314#chr4:55968027-55968093:-:333|0|chr4:55968023-55968093:+:549#chr4:55967974-55968185:+:321#chr4:55967975-55968186:-:307#chr4:55967939-55968092:+:437#chr4:55967938-55968091:-:413#chr4:55968027-55968093:-:436|0\n")
            variations.write("1133\t1130\t0.9973521624007061\tNC_000004.11\t1809120\tG\t0|1130|0|0\tG\tA/C/T\t1.0\t0|0|0|0\t365|225|279|261\t0|0|0|0\t0|0|0|0\t0|chr4:1809074-1809319:-:58#chr4:1809007-1809187:+:519#chr4:1809071-1809315:+:104#chr4:1809006-1809186:-:435|0|0\n")
            variations.write("1132\t1127\t0.9955830388692579\tNC_000004.11\t1809121\tT\t0|0|2|1127\tT\tC\t0.9982285208148804\t0|0|0|0\t0|0|0|0\t0|1|1|0\t362|224|280|261\t0|0|chr4:1809006-1809186:-:2|chr4:1809074-1809319:-:59#chr4:1809007-1809187:+:515#chr4:1809071-1809315:+:105#chr4:1809006-1809186:-:434\n")
            variations.write("1135\t1119\t0.9859030837004406\tNC_000004.11\t1809122\tG\t0|1119|0|0\tG\tA/C/T\t1.0\t0|0|0|0\t358|222|280|259\t0|0|0|0\t0|0|0|0\t0|chr4:1809074-1809319:-:59#chr4:1809007-1809187:+:512#chr4:1809071-1809315:+:102#chr4:1809006-1809186:-:432|0|0\n")
            variations.write("1131\t1114\t0.9849690539345711\tNC_000004.11\t1809123\tT\t0|0|1|1114\tT\tC\t0.9991031390134529\t0|0|0|0\t0|0|0|0\t1|0|0|0\t360|222|276|256\t0|0|chr4:1809071-1809315:+:1|chr4:1809074-1809319:-:59#chr4:1809007-1809187:+:513#chr4:1809071-1809315:+:100#chr4:1809006-1809186:-:428\n")
            variations.write("1133\t1016\t0.8967343336275375\tNC_000004.11\t1809124\tG\t0|1016|0|0\tG\tA/C/T\t1.0\t0|0|0|0\t321|204|258|233\t0|0|0|0\t0|0|0|0\t0|chr4:1809074-1809319:-:57#chr4:1809007-1809187:+:462#chr4:1809071-1809315:+:90#chr4:1809006-1809186:-:397|0|0\n")
            variations.write("1139\t1021\t0.8964003511852502\tNC_000004.11\t1809125\tT\t0|0|1|1021\tT\tC\t0.9990215264187867\t0|0|0|0\t0|0|0|0\t0|0|1|0\t321|205|260|235\t0|0|chr4:1809006-1809186:-:1|chr4:1809074-1809319:-:56#chr4:1809007-1809187:+:464#chr4:1809071-1809315:+:90#chr4:1809006-1809186:-:400\n")
            variations.write("1131\t1012\t0.8947833775419982\tNC_000004.11\t1809126\tG\t1|1012|0|0\tG\tA\t0.9990128331688055\t0|0|0|1\t319|205|258|230\t0|0|0|0\t0|0|0|0\tchr4:1809007-1809187:+:1|chr4:1809074-1809319:-:57#chr4:1809007-1809187:+:457#chr4:1809071-1809315:+:90#chr4:1809006-1809186:-:397|0|0\n")
            variations.write("1115\t991\t0.8887892376681614\tNC_000004.11\t1809127\tC\t1|0|991|5\tC\tT\t0.9949799196787149\t0|1|0|0\t0|0|0|0\t316|198|252|225\t1|2|0|2\t0|0|chr4:1809074-1809319:-:56#chr4:1809007-1809187:+:450#chr4:1809071-1809315:+:89#chr4:1809006-1809186:-:386|chr4:1809007-1809187:+:3#chr4:1809006-1809186:-:2\n")
            variations.write("1139\t126\t0.1106233538191396\tNC_000004.11\t1809128\tG\t0|126|0|0\tG\tA/C/T\t1.0\t0|0|0|0\t48|22|26|30\t0|0|0|0\t0|0|0|0\t0|chr4:1809074-1809319:-:2#chr4:1809007-1809187:+:61#chr4:1809071-1809315:+:16#chr4:1809006-1809186:-:43|0|0\n")
            variations.write("1138\t21\t0.01845342706502636\tNC_000004.11\t1809129\tT\t0|0|104|21\tC\tT\t0.832\t0|0|0|0\t0|0|0|0\t40|18|18|28\t8|4|6|3\t0|0|chr4:1809074-1809319:-:2#chr4:1809007-1809187:+:55#chr4:1809071-1809315:+:12#chr4:1809006-1809186:-:32|chr4:1809007-1809187:+:7#chr4:1809071-1809315:+:4#chr4:1809006-1809186:-:9\n")
            variations.write("293\t288\t0.9829351535836177\tNC_000002.11\t212543689\tC\t0|0|288|0\tC\tA/G/T\t1.0\t0|0|0|0\t0|0|0|0\t17|13|256|2\t0|0|0|0\t0|0|chr2:212543655-212543803:+:19#chr2:212543655-212543803:-:224|0\n")
            #variations.write("1818\t1765\t0.9708470847084708\tNC_000002.11\t29420596\tC\t0|0|1765|0\tC\tA/G/T\t1.0\t0|0|0|0\t0|0|0|0\t261|638|491|375\t0|0|0|0\t0|0|chr2:29420443-29420678:-:154#chr2:29420443-29420676:+:118#chr2:29420549-29420694:-:975#chr2:29420550-29420694:+:516|0\n")
            #variations.write("1822\t1766\t0.969264544456641\tNC_000002.11\t29420597\tA\t1766|3|0|0\tA\tG\t0.998304126625212\t262|639|489|376\t0|1|1|1\t0|0|0|0\t0|0|0|0\tchr2:29420443-29420678:-:154#chr2:29420443-29420676:+:115#chr2:29420549-29420694:-:974#chr2:29420550-29420694:+:521|chr2:29420549-29420694:-:2#chr2:29420550-29420694:+:1|0|0\n")
        with open(os.path.join(self.tempdir,"deletions"),'w') as deletions:
            deletions.write("243\tNC_000007.13\t140508151\t15\t0.06172839506172839\t1(-2,0)|4(0,2)|1(-1,0)|4(0,1)|1(-1,2)|1(0,6)|1(0,5)|2(0,3)\t0|15|0|0\t(-1,0):chr7:140508026-140508255:-:1|(0,5):chr7:140508026-140508255:-:1|(-2,0):chr7:140508026-140508255:-:1|(0,3):chr7:140508026-140508255:-:1|(0,1):chr7:140508026-140508255:-:4|(0,2):chr7:140508026-140508255:-:4|(-1,2):chr7:140508026-140508255:-:1|(0,6):chr7:140508026-140508255:-:1\n")
            deletions.write("91\tNC_000007.13\t140499713\t2\t0.02197802197802198\t2(0,0)\t0|1|1|0\t0\n")
            deletions.write("1115\tNC_000004.11\t1809127\t118\t0.10582959641255606\t13(-5,0)|2(-16,5)|101(-3,0)|1(-1,0)|1(-7,0)\t45|21|24|28\t(-3,0):chr4:1809074-1809319:-:2#(-3,0):chr4:1809007-1809187:+:52#(-3,0):chr4:1809071-1809315:+:12#(-3,0):chr4:1809006-1809186:-:32|(-16,5):chr4:1809006-1809186:-:2|(-1,0):chr4:1809006-1809186:-:1|(-7,0):chr4:1809071-1809315:+:1|(-5,0):chr4:1809007-1809187:+:4#(-5,0):chr4:1809071-1809315:+:3#(-5,0):chr4:1809006-1809186:-:6\n")
            deletions.write("1139\tNC_000004.11\t1809128\t1013\t0.8893766461808604\t2(-17,4)|917(0,1)|7(0,5)|87(0,3)\t320|204|259|230\t(0,5):chr4:1809007-1809187:+:4#(0,5):chr4:1809071-1809315:+:1#(0,5):chr4:1809006-1809186:-:2|(-17,4):chr4:1809006-1809186:-:2|(0,3):chr4:1809074-1809319:-:4#(0,3):chr4:1809007-1809187:+:38#(0,3):chr4:1809071-1809315:+:5#(0,3):chr4:1809006-1809186:-:40|(0,1):chr4:1809074-1809319:-:53#(0,1):chr4:1809007-1809187:+:416#(0,1):chr4:1809071-1809315:+:84#(0,1):chr4:1809006-1809186:-:354\n")
            deletions.write("1138\tNC_000004.11\t1809129\t1013\t0.8901581722319859\t2(-18,3)|917(-1,0)|87(-1,2)|7(-1,4)\t320|204|259|230\t(-18,3):chr4:1809006-1809186:-:2|(-1,0):chr4:1809074-1809319:-:53#(-1,0):chr4:1809007-1809187:+:416#(-1,0):chr4:1809071-1809315:+:84#(-1,0):chr4:1809006-1809186:-:354|(-1,2):chr4:1809074-1809319:-:4#(-1,2):chr4:1809007-1809187:+:38#(-1,2):chr4:1809071-1809315:+:5#(-1,2):chr4:1809006-1809186:-:40|(-1,4):chr4:1809007-1809187:+:4#(-1,4):chr4:1809071-1809315:+:1#(-1,4):chr4:1809006-1809186:-:2\n")
        with open(os.path.join(self.tempdir,"insertions"),'w') as insertions:
            insertions.write("353\tNC_000002.11\t212543689\t6\t0.0169971671388102\t4A|2AA\t0|5|0|1\tAA:chr2:212543655-212543803:-:2|A:chr2:212543655-212543803:+:1#A:chr2:212543655-212543803:-:3")

        self.nc_to_chr = {
            'NC_000001.10': '1',
            'NC_000002.11': '2',
            'NC_000003.11': '3',
            'NC_000004.11': '4',
            'NC_000007.13': '7',
            'NC_000010.10': '10',
            'NC_000012.11': '12',
            'NC_000014.8': '14',
            'NC_000015.9': '15',
            'NC_000017.10': '17',
            'NC_000019.9': '19'}
        self.ref = dict([(key,{}) for key in self.nc_to_chr])

    def tearDown(self):
        # delete fixtures
        shutil.rmtree(self.tempdir)
        pass

    def test_ref_import(self):

        from scripts.lib.common.data.parser.jsnpmania import extract_ref_variant_info
        from copy import deepcopy
        ref = deepcopy(self.ref)
        with open(os.path.join(self.tempdir, 'variations'),'r') as ref_lines:
            for line in ref_lines:
                ref = extract_ref_variant_info(ref,line)

        expected_ref = {
            'NC_000001.10': {
                '115252473': {'reference': 'T', 'depth': 76}},
            'NC_000003.11': {},
            'NC_000010.10': {
                '123256344': {'reference': 'T', 'depth': 1601},
                '123256345': {'reference': 'T', 'depth': 795},
                '123256346': {'reference': 'A', 'depth': 567},
                '123239112': {'reference': 'G', 'depth': 700}},
            'NC_000012.11': {},
            'NC_000014.8': {
                '105258892': {'reference': 'C', 'depth': 3532},
                '105258893': {'reference': 'A', 'depth': 1},
                '105258894': {'reference': 'T', 'depth': 6434}},
            'NC_000015.9': {},
            'NC_000017.10': {},
            'NC_000019.9': {},
            'NC_000002.11': {'212543689': {'reference': 'C', 'depth': 288}},
                             #'29420596': {'reference': 'C', 'depth': 1765},
                             #'29420597': {'reference': 'A', 'depth': 1766}},
            'NC_000004.11': {'1809120': {'reference': 'G', 'depth': 1130},
                             '1809121': {'reference': 'T', 'depth': 1127},
                             '1809122': {'reference': 'G', 'depth': 1119},
                             '1809123': {'reference': 'T', 'depth': 1114},
                             '1809124': {'reference': 'G', 'depth': 1016},
                             '1809125': {'reference': 'T', 'depth': 1021},
                             '1809126': {'reference': 'G', 'depth': 1012},
                             '1809127': {'reference': 'C', 'depth': 991},
                             '1809128': {'reference': 'G', 'depth': 126},
                             '1809129': {'reference': 'T', 'depth': 21},
                             '55968053': {'reference': 'A', 'depth': 1711}},
            'NC_000007.13': {'140499713': {'reference': 'A', 'depth': 89},
                             '140508149': {'reference': 'T', 'depth': 216},
                             '140508150': {'reference': 'C', 'depth': 225},
                             '140508151': {'reference': 'C', 'depth': 221},
                             '140508152': {'reference': 'A', 'depth': 205},
                             '140508153': {'reference': 'A', 'depth': 217},
                             '140508154': {'reference': 'A', 'depth': 238},
                             '140508155': {'reference': 'A', 'depth': 259},
                             '140508156': {'reference': 'A', 'depth': 269},
                             '140508157': {'reference': 'A', 'depth': 279}}}
        for chr_key, values in ref.items():
            for position in ref[chr_key]:
                for key, value in  ref[chr_key][position].items():
                    self.assertEqual(expected_ref[chr_key][position][key], value)

        ref = deepcopy(self.ref)
        with open(os.path.join(self.tempdir, 'variations'),'r') as ref_lines:
            for line in ref_lines:
                ref = extract_ref_variant_info(ref,line,1)

        expected_ref = {
            'NC_000001.10': {
                '115252473': {'reference': 'T', 'depth': 76,
                    'amp+': 1, 'amp-': 1, 'ampInfo': 'chr1:115252332-115252647:+:43#chr1:115252332-115252648:-:33'}},
            'NC_000003.11': {},
            'NC_000010.10': {
                '123256344': {'reference': 'T', 'depth': 1601,
                    'amp+': 1, 'amp-': 1, 'ampInfo': 'chr10:123256234-123256344:+:416#chr10:123256234-123256346:-:588'},
                '123256345': {'reference': 'T', 'depth': 795,
                    'amp+': 1, 'amp-': 1, 'ampInfo': 'chr10:123256234-123256344:+:1#chr10:123256234-123256346:-:587'},
                '123256346': {'reference': 'A', 'depth': 567,
                    'amp+': 0, 'amp-': 1, 'ampInfo': 'chr10:123256234-123256346:-:564'},
                '123239112': {'reference': 'G', 'depth': 700,
                    'amp+': 2, 'amp-': 2, 'ampInfo': 'chr10:123239059-123239249:-:40#chr10:123239055-123239245:+:406#chr10:123238910-123239122:+:192#chr10:123238910-123239121:-:47'}
            },
            'NC_000012.11': {},
            'NC_000014.8': {
                '105258892': {'reference': 'C', 'depth': 3532,
                    'amp+': 3, 'amp-': 3, 'ampInfo': 'chr14:105258859-105258986:+:1544#chr14:105258855-105258982:-:1010#chr14:105258810-105258944:-:324#chr14:105258806-105258944:+:363#chr14:105258778-105259019:-:52#chr14:105258774-105259016:+:228'},
                '105258893': {'reference': 'A', 'depth': 1,
                    'amp+': 0, 'amp-': 1, 'ampInfo': 'chr14:105258855-105258982:-:1'},
                '105258894': {'reference': 'T', 'depth': 6434,
                    'amp+': 3, 'amp-': 3, 'ampInfo': 'chr14:105258859-105258986:+:2731#chr14:105258855-105258982:-:1917#chr14:105258810-105258944:-:595#chr14:105258806-105258944:+:689#chr14:105258778-105259019:-:96#chr14:105258774-105259016:+:368'}
            },
            'NC_000015.9': {},
            'NC_000017.10': {},
            'NC_000019.9': {},
            'NC_000002.11': {'212543689': {'reference': 'C', 'depth': 288,
                                'amp+': 1, 'amp-': 1, 'ampInfo': 'chr2:212543655-212543803:+:19#chr2:212543655-212543803:-:224'}},
                             #'29420596': {'reference': 'C', 'depth': 1765,
                             #   'amp+': 2, 'amp-': 2, 'ampInfo': 'chr2:29420443-29420678:-:154#chr2:29420443-29420676:+:118#chr2:29420549-29420694:-:975#chr2:29420550-29420694:+:516'},
                             #'29420597': {'reference': 'A', 'depth': 1766,
                             #   'amp+': 2, 'amp-': 2, 'ampInfo': 'chr2:29420443-29420678:-:154#chr2:29420443-29420676:+:115#chr2:29420549-29420694:-:974#chr2:29420550-29420694:+:521'}},
            'NC_000004.11': {'1809120': {'reference': 'G', 'depth': 1130,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:58#chr4:1809007-1809187:+:519#chr4:1809071-1809315:+:104#chr4:1809006-1809186:-:435'},
                             '1809121': {'reference': 'T', 'depth': 1127,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:59#chr4:1809007-1809187:+:515#chr4:1809071-1809315:+:105#chr4:1809006-1809186:-:434'},
                             '1809121': {'reference': 'T', 'depth': 1127,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:59#chr4:1809007-1809187:+:515#chr4:1809071-1809315:+:105#chr4:1809006-1809186:-:434'},
                             '1809122': {'reference': 'G', 'depth': 1119,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:59#chr4:1809007-1809187:+:512#chr4:1809071-1809315:+:102#chr4:1809006-1809186:-:432'},
                             '1809123': {'reference': 'T', 'depth': 1114,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:59#chr4:1809007-1809187:+:513#chr4:1809071-1809315:+:100#chr4:1809006-1809186:-:428'},
                             '1809124': {'reference': 'G', 'depth': 1016,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:57#chr4:1809007-1809187:+:462#chr4:1809071-1809315:+:90#chr4:1809006-1809186:-:397'},
                             '1809125': {'reference': 'T', 'depth': 1021,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:56#chr4:1809007-1809187:+:464#chr4:1809071-1809315:+:90#chr4:1809006-1809186:-:400'},
                             '1809126': {'reference': 'G', 'depth': 1012,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:57#chr4:1809007-1809187:+:457#chr4:1809071-1809315:+:90#chr4:1809006-1809186:-:397'},
                             '1809127': {'reference': 'C', 'depth': 991,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:56#chr4:1809007-1809187:+:450#chr4:1809071-1809315:+:89#chr4:1809006-1809186:-:386'},
                             '1809128': {'reference': 'G', 'depth': 126,
                                'amp+': 2, 'amp-': 2, 'ampInfo': 'chr4:1809074-1809319:-:2#chr4:1809007-1809187:+:61#chr4:1809071-1809315:+:16#chr4:1809006-1809186:-:43'},
                             '1809129': {'reference': 'T', 'depth': 21,
                                'amp+': 2, 'amp-': 1, 'ampInfo': 'chr4:1809007-1809187:+:7#chr4:1809071-1809315:+:4#chr4:1809006-1809186:-:9'},
                            '55968053': {'reference': 'A', 'depth': 1711,
                                'amp+': 3, 'amp-': 3, 'ampInfo': 'chr4:55968023-55968093:+:384#chr4:55967974-55968185:+:243#chr4:55967975-55968186:-:188#chr4:55967939-55968092:+:232#chr4:55967938-55968091:-:314#chr4:55968027-55968093:-:333'}},
            'NC_000007.13': {'140499713': {'reference': 'A', 'depth': 89,
                                'amp+': 1, 'amp-': 1, 'ampInfo': 'chr7:140499687-140500177:-:33#chr7:140499687-140500177:+:44'},
                             '140508149': {'reference': 'T', 'depth': 216,
                                'amp+': 1, 'amp-': 2, 'ampInfo': 'chr7:140508000-140508412:+:36#chr7:140508001-140508412:-:19#chr7:140508026-140508255:-:155'},
                             '140508150': {'reference': 'C', 'depth': 225,
                                'amp+': 1, 'amp-': 2, 'ampInfo': 'chr7:140508000-140508412:+:34#chr7:140508001-140508412:-:18#chr7:140508026-140508255:-:166'},
                             '140508151': {'reference': 'C', 'depth': 221,
                                'amp+': 1, 'amp-': 2, 'ampInfo': 'chr7:140508000-140508412:+:26#chr7:140508001-140508412:-:17#chr7:140508026-140508255:-:172'},
                             '140508152': {'reference': 'A', 'depth': 205,
                                'amp+': 0, 'amp-': 2, 'ampInfo': 'chr7:140508001-140508412:-:19#chr7:140508026-140508255:-:177'},
                             '140508153': {'reference': 'A', 'depth': 217,
                                'amp+': 0, 'amp-': 1, 'ampInfo': 'chr7:140508026-140508255:-:206'},
                             '140508154': {'reference': 'A', 'depth': 238,
                                'amp+': 0, 'amp-': 1, 'ampInfo': 'chr7:140508026-140508255:-:227'},
                             '140508155': {'reference': 'A', 'depth': 259,
                                'amp+': 0, 'amp-': 1, 'ampInfo': 'chr7:140508026-140508255:-:247'},
                             '140508156': {'reference': 'A', 'depth': 269,
                                'amp+': 0, 'amp-': 1, 'ampInfo': 'chr7:140508026-140508255:-:258'},
                             '140508157': {'reference': 'A', 'depth': 279,
                                'amp+': 0, 'amp-': 1, 'ampInfo': 'chr7:140508026-140508255:-:267'}}}

        for chr_key in ref:
            self.assertEqual(len(expected_ref[chr_key].keys()), len(ref[chr_key]))
            for position in ref[chr_key]:
                for key, value in  ref[chr_key][position].items():
                    self.assertEqual(expected_ref[chr_key][position][key], value)

        self.assertEqual(extract_ref_variant_info({'NC_000007.13': {}},"36\t34\t0.9444444444444444\tNC_000007.13\t140487214\tG\t34|2|0|0\tG\tT\t0.9444444444444444\t0|8|26|0\t0|1|1|0\t0|0|0|0\t0|0|0|0\tchr7:140487208-140487699:-:18|0|0|0",1),
                         {'NC_000007.13': {'140487214': {'reference': 'G', 'depth': 2, 'amp+': 0, 'amp-': 0, 'ampInfo': '0'}}})
        self.assertRaises(AttributeError,
                            extract_ref_variant_info, {},"36\t34\t0.9444444444444444\tNC_000007.13\t140487214\tG\t34|2|0|0\tA\tT\t0.9444444444444444\t0|8|26|0\t0|1|1|0\t0|0|0|0\t0|0|0|0\tchr7:140487208-140487699:-:18|0|0|0",1)

    def test_get_major_vaf(self):
        from scripts.lib.common.data.parser.jsnpmania import _get_major_vaf
        from collections import OrderedDict
        self.assertEqual(('NC_000003.11#178919375#178919375#C#A', 0.169811320754717),
            _get_major_vaf(
            OrderedDict([(
                'NC_000003.11#178919375#178919375#C#T',
                      '3\t178919375\t178919375\tC\tT\tcomments: sample=17-1891 variantAlleleRatio=0.0188679245283019 alleleFreq=42,1 readDepth=53 Tumor_A=0|3|0|6 Tumor_G=0|0|0|0 Tumor_C=0|21|0|21 Tumor_T=0|0|0|1 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=chr3:178919054-178919447:+:1 Tumor_ref_ampliconInfo=chr3:178919156-178919524:+:1#chr3:178919054-178919447:-:15#chr3:178919054-178919447:+:20#chr3:178919156-178919524:-:6'),
                 ('NC_000003.11#178919375#178919375#C#A',
                      '3\t178919375\t178919375\tC\tA\tcomments: sample=17-1891 variantAlleleRatio=0.169811320754717 alleleFreq=42,9 readDepth=53 Tumor_A=0|3|0|6 Tumor_G=0|0|0|0 Tumor_C=0|21|0|21 Tumor_T=0|0|0|1 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=chr3:178919156-178919524:+:3#chr3:178919054-178919447:-:2#chr3:178919054-178919447:+:3#chr3:178919156-178919524:-:1 Tumor_ref_ampliconInfo=chr3:178919156-178919524:+:1#chr3:178919054-178919447:-:15#chr3:178919054-178919447:+:20#chr3:178919156-178919524:-:6')])))


    def test_variants_import(self):
        from scripts.lib.common.data.parser.jsnpmania import extract_ref_variant_info, extract_major_snv_allele, _column_converter_snv
        from copy import deepcopy

        ref = deepcopy(self.ref)
        with open(os.path.join(self.tempdir, 'variations'),'r') as ref_lines:
            for line in ref_lines:
                ref = extract_ref_variant_info(ref,line,20)

        variant = "4199\t1711\t0.40747797094546323\tNC_000004.11\t55968053" + \
                  "\tA\t1711|0|2487|0\tC\tA\t0.5924249642686994\t428|426|426|431" + \
                  "\t0|0|0|0\t654|590|575|668\t0|0|0|0" + \
                  "\tchr4:55968023-55968093:+:384#chr4:55967974-55968185:+:243#chr4:55967975-55968186:-:188#chr4:55967939-55968092:+:232#chr4:55967938-55968091:-:314#chr4:55968027-55968093:-:333|0|chr4:55968023-55968093:+:549#chr4:55967974-55968185:+:321#chr4:55967975-55968186:-:307#chr4:55967939-55968092:+:437#chr4:55967938-55968091:-:413#chr4:55968027-55968093:-:436|0"

        def allele_filter(allele_depth, columns): return int(columns[0]) >= 20 and int(allele_depth)/int(columns[0]) >= 0.01 and (1.0 - float(columns[_column_converter_snv['allele_ratio']])) >= 0.01
        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 0, allele_filter)
        self.maxDiff = None
        self.assertIn('NC_000004.11#55968053#55968053#A#C', info.keys())
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("4\t55968053\t55968053\tA\tC\tcomments: sample=sample1 variantAlleleRatio=0.5922838771135984 alleleFreq=1711,2487 readDepth=4199 Tumor_A=428|426|426|431 Tumor_G=0|0|0|0 Tumor_C=654|590|575|668 Tumor_T=0|0|0|0",
            info['NC_000004.11#55968053#55968053#A#C'])

        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 1,allele_filter)
        self.assertIn('NC_000004.11#55968053#55968053#A#C', info.keys())
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("4\t55968053\t55968053\tA\tC\tcomments: sample=sample1 variantAlleleRatio=0.5922838771135984 alleleFreq=1711,2487 readDepth=4199 Tumor_A=428|426|426|431 Tumor_G=0|0|0|0 Tumor_C=654|590|575|668 Tumor_T=0|0|0|0 Tumor_var_plusAmplicons=3 Tumor_var_minusAmplicons=3 Tumor_ref_plusAmplicons=3 Tumor_ref_minusAmplicons=3 Tumor_var_ampliconInfo=chr4:55968023-55968093:+:549#chr4:55967974-55968185:+:321#chr4:55967975-55968186:-:307#chr4:55967939-55968092:+:437#chr4:55967938-55968091:-:413#chr4:55968027-55968093:-:436 Tumor_ref_ampliconInfo=chr4:55968023-55968093:+:384#chr4:55967974-55968185:+:243#chr4:55967975-55968186:-:188#chr4:55967939-55968092:+:232#chr4:55967938-55968091:-:314#chr4:55968027-55968093:-:333", info['NC_000004.11#55968053#55968053#A#C'])

        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 5,allele_filter)
        self.assertIn('NC_000004.11#55968053#55968053#A#C', info.keys())
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("4\t55968053\t55968053\tA\tC\tcomments: sample=sample1 variantAlleleRatio=0.5922838771135984 alleleFreq=1711,2487 readDepth=4199 Tumor_A=428|426|426|431 Tumor_G=0|0|0|0 Tumor_C=654|590|575|668 Tumor_T=0|0|0|0 Tumor_var_plusAmplicons=3 Tumor_var_minusAmplicons=3 Tumor_ref_plusAmplicons=3 Tumor_ref_minusAmplicons=3 Tumor_var_ampliconInfo=chr4:55968023-55968093:+:549#chr4:55967974-55968185:+:321#chr4:55967975-55968186:-:307#chr4:55967939-55968092:+:437#chr4:55967938-55968091:-:413#chr4:55968027-55968093:-:436 Tumor_ref_ampliconInfo=chr4:55968023-55968093:+:384#chr4:55967974-55968185:+:243#chr4:55967975-55968186:-:188#chr4:55967939-55968092:+:232#chr4:55967938-55968091:-:314#chr4:55968027-55968093:-:333", info['NC_000004.11#55968053#55968053#A#C'])

        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 380,allele_filter)
        self.assertIn('NC_000004.11#55968053#55968053#A#C', info.keys())
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("4\t55968053\t55968053\tA\tC\tcomments: sample=sample1 variantAlleleRatio=0.5922838771135984 alleleFreq=1711,2487 readDepth=4199 Tumor_A=428|426|426|431 Tumor_G=0|0|0|0 Tumor_C=654|590|575|668 Tumor_T=0|0|0|0 Tumor_var_plusAmplicons=2 Tumor_var_minusAmplicons=2 Tumor_ref_plusAmplicons=3 Tumor_ref_minusAmplicons=3 Tumor_var_ampliconInfo=chr4:55968023-55968093:+:549#chr4:55967974-55968185:+:321#chr4:55967975-55968186:-:307#chr4:55967939-55968092:+:437#chr4:55967938-55968091:-:413#chr4:55968027-55968093:-:436 Tumor_ref_ampliconInfo=chr4:55968023-55968093:+:384#chr4:55967974-55968185:+:243#chr4:55967975-55968186:-:188#chr4:55967939-55968092:+:232#chr4:55967938-55968091:-:314#chr4:55968027-55968093:-:333", info['NC_000004.11#55968053#55968053#A#C'])

        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 380, lambda variant_depth, columns: int(variant_depth) > 2487)
        self.assertEqual(0, len(info.keys()))
        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 380, lambda variant_depth, columns: int(variant_depth) >= 2487)
        self.assertEqual(1, len(info.keys()))

        variant = "78\t76\t0.9743589743589743\tNC_000001.10\t115252473\tT\t1|0|1|76\tT\tA/C\t0.987012987012987\t0|0|1|0\t0|0|0|0\t1|0|0|0\t43|0|33|0\tchr1:115252332-115252648:-:1|0|chr1:115252332-115252647:+:1|chr1:115252332-115252647:+:43#chr1:115252332-115252648:-:33"
        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 20, lambda variant_depth, columns: int(variant_depth)/int(columns[_column_converter_snv['depth']]) >= 0.01)
        self.assertEqual(2, len(info.keys()))
        self.assertEqual("1\t115252473\t115252473\tT\tC\tcomments: sample=sample1 variantAlleleRatio=0.0128205128205128 alleleFreq=76,1 readDepth=78 Tumor_A=0|0|1|0 Tumor_G=0|0|0|0 Tumor_C=1|0|0|0 Tumor_T=43|0|33|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=chr1:115252332-115252647:+:1 Tumor_ref_ampliconInfo=chr1:115252332-115252647:+:43#chr1:115252332-115252648:-:33", info['NC_000001.10#115252473#115252473#T#C'])
        self.assertEqual("1\t115252473\t115252473\tT\tA\tcomments: sample=sample1 variantAlleleRatio=0.0128205128205128 alleleFreq=76,1 readDepth=78 Tumor_A=0|0|1|0 Tumor_G=0|0|0|0 Tumor_C=1|0|0|0 Tumor_T=43|0|33|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=chr1:115252332-115252648:-:1 Tumor_ref_ampliconInfo=chr1:115252332-115252647:+:43#chr1:115252332-115252648:-:33", info['NC_000001.10#115252473#115252473#T#A'])

        variant = "6463\t3532\t0.5464954355562432\tNC_000014.8\t105258892\tC\t4|0|3532|2927\tC\tT\t0.5468338752128813\t0|2|2|0\t0|0|0|0\t1077|691|695|1069\t842|627|629|829\tchr14:105258855-105258982:-:2#chr14:105258810-105258944:-:2|0|chr14:105258859-105258986:+:1544#chr14:105258855-105258982:-:1010#chr14:105258810-105258944:-:324#chr14:105258806-105258944:+:363#chr14:105258778-105259019:-:52#chr14:105258774-105259016:+:228|chr14:105258859-105258986:+:1193#chr14:105258855-105258982:-:916#chr14:105258810-105258944:-:274#chr14:105258806-105258944:+:323#chr14:105258778-105259019:-:47#chr14:105258774-105259016:+:145\n"
        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 20, allele_filter)
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("14\t105258892\t105258892\tC\tT\tcomments: sample=sample1 variantAlleleRatio=0.4528856568157202 alleleFreq=3532,2927 readDepth=6463 Tumor_A=0|2|2|0 Tumor_G=0|0|0|0 Tumor_C=1077|691|695|1069 Tumor_T=842|627|629|829 Tumor_var_plusAmplicons=3 Tumor_var_minusAmplicons=3 Tumor_ref_plusAmplicons=3 Tumor_ref_minusAmplicons=3 Tumor_var_ampliconInfo=chr14:105258859-105258986:+:1193#chr14:105258855-105258982:-:916#chr14:105258810-105258944:-:274#chr14:105258806-105258944:+:323#chr14:105258778-105259019:-:47#chr14:105258774-105259016:+:145 Tumor_ref_ampliconInfo=chr14:105258859-105258986:+:1544#chr14:105258855-105258982:-:1010#chr14:105258810-105258944:-:324#chr14:105258806-105258944:+:363#chr14:105258778-105259019:-:52#chr14:105258774-105259016:+:228",info["NC_000014.8#105258892#105258892#C#T"])
        variant = "6455\t1\t1.5491866769945777E-4\tNC_000014.8\t105258893\tA\t1|6454|0|0\tG\tA\t0.9998450813323005\t0|0|1|0\t1915|1318|1320|1901\t0|0|0|0\t0|0|0|0\tchr14:105258855-105258982:-:1|chr14:105258859-105258986:+:2737#chr14:105258855-105258982:-:1929#chr14:105258810-105258944:-:595#chr14:105258806-105258944:+:686#chr14:105258778-105259019:-:95#chr14:105258774-105259016:+:374|0|0\n"
        info = extract_major_snv_allele('sample1', ref, variant, self.nc_to_chr, 20, allele_filter)
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("14\t105258893\t105258893\tA\tG\tcomments: sample=sample1 variantAlleleRatio=0.9998450813323005 alleleFreq=1,6454 readDepth=6455 Tumor_A=0|0|1|0 Tumor_G=1915|1318|1320|1901 Tumor_C=0|0|0|0 Tumor_T=0|0|0|0 Tumor_var_plusAmplicons=3 Tumor_var_minusAmplicons=3 Tumor_ref_plusAmplicons=0 Tumor_ref_minusAmplicons=0 Tumor_var_ampliconInfo=chr14:105258859-105258986:+:2737#chr14:105258855-105258982:-:1929#chr14:105258810-105258944:-:595#chr14:105258806-105258944:+:686#chr14:105258778-105259019:-:95#chr14:105258774-105259016:+:374 Tumor_ref_ampliconInfo=chr14:105258855-105258982:-:1",info["NC_000014.8#105258893#105258893#A#G"])


    def test_deletion_import(self):
        from scripts.lib.common.data.parser.jsnpmania import extract_ref_variant_info, extract_deletions, _column_converter_indel, _update_deletion_data
        from copy import deepcopy

        ref = deepcopy(self.ref)
        with open(os.path.join(self.tempdir, 'variations'),'r') as ref_lines:
            for line in ref_lines:
                ref = extract_ref_variant_info(ref,line,1)

        self.maxDiff = None

        def deletion_filter(variant_depth, columns): return int(columns[_column_converter_indel['depth']]) >= 20 and int(variant_depth)/int(columns[_column_converter_indel['depth']]) >= 0.01
        deletion = "243\tNC_000007.13\t140508151\t15\t0.06172839506172839\t1(-2,0)|4(0,2)|1(-1,0)|4(0,1)|1(-1,2)|1(0,6)|1(0,5)|2(0,3)\t0|15|0|0\t(-1,0):chr7:140508026-140508255:-:1|(0,5):chr7:140508026-140508255:-:1|(-2,0):chr7:140508026-140508255:-:1|(0,3):chr7:140508026-140508255:-:1|(0,1):chr7:140508026-140508255:-:4|(0,2):chr7:140508026-140508255:-:4|(-1,2):chr7:140508026-140508255:-:1|(0,6):chr7:140508026-140508255:-:1"
        info = extract_deletions("sample1", ref, deletion, self.nc_to_chr, 0, deletion_filter)
        self.assertIn('NC_000007.13#140508151#140508153##-', info.keys())
        self.assertIn('NC_000007.13#140508151#140508152##-', info.keys())
        self.assertEqual(2, len(info.keys()))
        self.assertEqual("7\t140508151\t140508153\tCAA\t-\tcomments: sample=sample1 variantAlleleRatio=0.01646090534979424 alleleFreq=221,4 readDepth=243 Tumor_Del=0|15|0|0",info['NC_000007.13#140508151#140508153##-'])
        self.assertEqual("7\t140508151\t140508152\tCA\t-\tcomments: sample=sample1 variantAlleleRatio=0.01646090534979424 alleleFreq=221,4 readDepth=243 Tumor_Del=0|15|0|0",info['NC_000007.13#140508151#140508152##-'])
        info = extract_deletions("sample1", ref, deletion, self.nc_to_chr, 5, deletion_filter)
        self.assertEqual("7\t140508151\t140508153\tCAA\t-\tcomments: sample=sample1 variantAlleleRatio=0.01646090534979424 alleleFreq=221,4 readDepth=243 Tumor_Del=0|15|0|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=(0,2):chr7:140508026-140508255:-:4 Tumor_ref_ampliconInfo=chr7:140508000-140508412:+:26#chr7:140508001-140508412:-:17#chr7:140508026-140508255:-:172",info['NC_000007.13#140508151#140508153##-'])
        self.assertEqual("7\t140508151\t140508152\tCA\t-\tcomments: sample=sample1 variantAlleleRatio=0.01646090534979424 alleleFreq=221,4 readDepth=243 Tumor_Del=0|15|0|0 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=(0,1):chr7:140508026-140508255:-:4 Tumor_ref_ampliconInfo=chr7:140508000-140508412:+:26#chr7:140508001-140508412:-:17#chr7:140508026-140508255:-:172",info['NC_000007.13#140508151#140508152##-'])

        def deletion_filter(variant_depth, columns): return int(columns[_column_converter_indel['depth']]) >= 20 and int(variant_depth)/int(columns[_column_converter_indel['depth']]) >= 0.02
        info = extract_deletions("sample1", ref, deletion, self.nc_to_chr, 5, deletion_filter)
        self.assertEqual(0, len(info.keys()))
        deletion = "91\tNC_000007.13\t140499713\t2\t0.02197802197802198\t2(0,0)\t0|1|1|0\t0"
        info = extract_deletions("sample1", ref, deletion, self.nc_to_chr, 5, deletion_filter)
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("7\t140499713\t140499713\tA\t-\tcomments: sample=sample1 variantAlleleRatio=0.02197802197802198 alleleFreq=89,2 readDepth=91 Tumor_Del=0|1|1|0 Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=- Tumor_ref_ampliconInfo=chr7:140499687-140500177:-:33#chr7:140499687-140500177:+:44",info['NC_000007.13#140499713#140499713##-'])

        deletion_2 = "1135\tNC_000004.11\t1809122\t16\t0.014096916299559472\t1(-2,5)|13(0,5)\t6|4|4|2\t(-2,5):chr4:1809071-1809315:+:1|(0,5):chr4:1809007-1809187:+:4#(0,5):chr4:1809071-1809315:+:3#(0,5):chr4:1809006-1809186:-:6"
        info = extract_deletions("sample1", ref, deletion_2, self.nc_to_chr, 1, None)
        self.assertEqual(
            "4\t1809122\t1809127\tGTGTGC\t-\tcomments: sample=sample1 variantAlleleRatio=0.01145374449339207 alleleFreq=1119,13 readDepth=1135 Tumor_Del=6|4|4|2 Tumor_var_plusAmplicons=2 Tumor_var_minusAmplicons=1 Tumor_ref_plusAmplicons=2 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=(0,5):chr4:1809007-1809187:+:4#(0,5):chr4:1809071-1809315:+:3#(0,5):chr4:1809006-1809186:-:6 Tumor_ref_ampliconInfo=chr4:1809074-1809319:-:59#chr4:1809007-1809187:+:512#chr4:1809071-1809315:+:102#chr4:1809006-1809186:-:432",
            info['NC_000004.11#1809122#1809127##-'])
        deletion_3 = "1133\tNC_000004.11\t1809124\t117\t0.10326566637246248\t13(-2,3)|1(-4,3)|101(0,3)\t45|21|23|28\t(-2,3):chr4:1809007-1809187:+:4#(-2,3):chr4:1809071-1809315:+:3#(-2,3):chr4:1809006-1809186:-:6|(0,3):chr4:1809074-1809319:-:2#(0,3):chr4:1809007-1809187:+:52#(0,3):chr4:1809071-1809315:+:12#(0,3):chr4:1809006-1809186:-:32|(-4,3):chr4:1809071-1809315:+:1"
        info = _update_deletion_data(info, extract_deletions("sample1", ref, deletion_3, self.nc_to_chr, 1, None))
        self.assertEqual(
            "4\t1809122\t1809127\tGTGTGC\t-\tcomments: sample=sample1 variantAlleleRatio=0.01147396293027361 alleleFreq=1016,13 readDepth=1133 Tumor_Del=45|21|23|28 Tumor_var_plusAmplicons=2 Tumor_var_minusAmplicons=1 Tumor_ref_plusAmplicons=2 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=(-2,3):chr4:1809007-1809187:+:4#(-2,3):chr4:1809071-1809315:+:3#(-2,3):chr4:1809006-1809186:-:6 Tumor_ref_ampliconInfo=chr4:1809074-1809319:-:57#chr4:1809007-1809187:+:462#chr4:1809071-1809315:+:90#chr4:1809006-1809186:-:397",
            info['NC_000004.11#1809122#1809127##-'])
        deletion_4 = "1115\tNC_000004.11\t1809127\t118\t0.10582959641255606\t13(-5,0)|101(-3,0)|1(-1,0)|1(-7,0)\t45|21|24|28\t(-3,0):chr4:1809074-1809319:-:2#(-3,0):chr4:1809007-1809187:+:52#(-3,0):chr4:1809071-1809315:+:12#(-3,0):chr4:1809006-1809186:-:32|(-1,0):chr4:1809006-1809186:-:1|(-7,0):chr4:1809071-1809315:+:1|(-5,0):chr4:1809007-1809187:+:4#(-5,0):chr4:1809071-1809315:+:3#(-5,0):chr4:1809006-1809186:-:6\n"
        info = _update_deletion_data(info, extract_deletions("sample1", ref, deletion_4, self.nc_to_chr, 5, None))
        self.assertEqual("4\t1809122\t1809127\tGTGTGC\t-\tcomments: sample=sample1 variantAlleleRatio=0.011659192825112108 alleleFreq=991,13 readDepth=1115 Tumor_Del=45|21|24|28 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=1 Tumor_ref_plusAmplicons=2 Tumor_ref_minusAmplicons=2 Tumor_var_ampliconInfo=(-5,0):chr4:1809007-1809187:+:4#(-5,0):chr4:1809071-1809315:+:3#(-5,0):chr4:1809006-1809186:-:6 Tumor_ref_ampliconInfo=chr4:1809074-1809319:-:56#chr4:1809007-1809187:+:450#chr4:1809071-1809315:+:89#chr4:1809006-1809186:-:386",
                         info['NC_000004.11#1809122#1809127##-'])

    def test_count_amplicons(self):
        from scripts.lib.common.data.parser.jsnpmania import _count_amplicons
        self.assertEqual(_count_amplicons("chr1:162722740-162722971:+:1#chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195", 2),(1, 1))
        self.assertEqual(_count_amplicons("chr1:162722740-162722971:+:1#chr1:162722740-162722971:-:241#chr1:162722740-162722971:+:195"),(2, 1))
        self.assertEqual(_count_amplicons("(0,0):chr1:162722740-162722971:+:1"), (1, 0))

    def test_count_amplicons(self):
        from scripts.lib.common.data.parser.jsnpmania import _get_ref_amplicon_information
        self.assertEqual(_get_ref_amplicon_information(['','chr7:140534338-140534738:+:14#chr7:140534619-140534946:-:127#chr7:140534531-140534841:-:99#chr7:140534532-140534841:+:81|0|0|0'], 0, 1),('chr7:140534338-140534738:+:14#chr7:140534619-140534946:-:127#chr7:140534531-140534841:-:99#chr7:140534532-140534841:+:81', 2, 2))
        self.assertEqual(_get_ref_amplicon_information(['','0|0|0|0'],0,1),('0',0,0))

    def test_insertion_import(self):
        from scripts.lib.common.data.parser.jsnpmania import extract_ref_variant_info, extract_insertion, _column_converter_indel
        from copy import deepcopy

        ref = deepcopy(self.ref)
        with open(os.path.join(self.tempdir, 'variations'),'r') as ref_lines:
            for line in ref_lines:
                ref = extract_ref_variant_info(ref,line,1)

        self.maxDiff = None

        def ifilter(variant_depth, columns):  return int(columns[_column_converter_indel['depth']]) >= 20 and int(variant_depth)/int(columns[_column_converter_indel['depth']]) >= 0.01
        insert = "353\tNC_000002.11\t212543689\t6\t0.0169971671388102\t4A|2AA\t0|5|0|1\tAA:chr2:212543655-212543803:-:2|A:chr2:212543655-212543803:+:1#A:chr2:212543655-212543803:-:3"
        info = extract_insertion('sample1', ref, insert, self.nc_to_chr, insertion_filter=ifilter)
        self.assertIn('NC_000002.11#212543689#212543689#-#A', info.keys())
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("2\t212543689\t212543689\t-\tA\tcomments: sample=sample1 variantAlleleRatio=0.0113314447592068 alleleFreq=288,4 readDepth=353 Tumor_Ins=0|5|0|1",info['NC_000002.11#212543689#212543689#-#A'])
        info = extract_insertion('sample1', ref, insert, self.nc_to_chr, 1, ifilter)
        self.assertIn('NC_000002.11#212543689#212543689#-#A', info.keys())
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("2\t212543689\t212543689\t-\tA\tcomments: sample=sample1 variantAlleleRatio=0.0113314447592068 alleleFreq=288,4 readDepth=353 Tumor_Ins=0|5|0|1 Tumor_var_plusAmplicons=1 Tumor_var_minusAmplicons=1 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=A:chr2:212543655-212543803:+:1#A:chr2:212543655-212543803:-:3 Tumor_ref_ampliconInfo=chr2:212543655-212543803:+:19#chr2:212543655-212543803:-:224",info['NC_000002.11#212543689#212543689#-#A'])
        info = extract_insertion('sample1', ref, insert, self.nc_to_chr, 5, ifilter)
        self.assertIn('NC_000002.11#212543689#212543689#-#A', info.keys())
        self.assertEqual(1, len(info.keys()))
        self.assertEqual("2\t212543689\t212543689\t-\tA\tcomments: sample=sample1 variantAlleleRatio=0.0113314447592068 alleleFreq=288,4 readDepth=353 Tumor_Ins=0|5|0|1 Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=A:chr2:212543655-212543803:+:1#A:chr2:212543655-212543803:-:3 Tumor_ref_ampliconInfo=chr2:212543655-212543803:+:19#chr2:212543655-212543803:-:224",info['NC_000002.11#212543689#212543689#-#A'])
        def ifilter( variant_depth, columns): return int(columns[_column_converter_indel['depth']]) >= 20 and int(variant_depth)/int(columns[_column_converter_indel['depth']]) >= 0.02
        info = extract_insertion('sample1', ref, insert, self.nc_to_chr, 5, ifilter)
        self.assertEqual(0, len(info.keys()))

    def test_multi_bp_import(self):
        from scripts.lib.common.data.parser.jsnpmania import extract_ref_variant_info, extract_major_snv_allele, _column_converter_snv, generate_multi_bp_variants
        from copy import deepcopy
        ref = deepcopy(self.ref)
        with open(os.path.join(self.tempdir, 'variations'),'r') as ref_lines:
            for line in ref_lines:
                ref = extract_ref_variant_info(ref,line,20)
        variants = dict()

        def allele_filter(allele_depth, columns): return int(columns[0]) >= 20 and int(allele_depth)/int(columns[0]) >= 0.01 and (1.0 - float(columns[_column_converter_snv['allele_ratio']])) >= 0.01

        variant_1 = "6463\t3532\t0.5464954355562432\tNC_000014.8\t105258892\tC\t4|0|3532|2927\tC\tT\t0.5468338752128813\t0|2|2|0\t0|0|0|0\t1077|691|695|1069\t842|627|629|829\tchr14:105258855-105258982:-:2#chr14:105258810-105258944:-:2|0|chr14:105258859-105258986:+:1544#chr14:105258855-105258982:-:1010#chr14:105258810-105258944:-:324#chr14:105258806-105258944:+:363#chr14:105258778-105259019:-:52#chr14:105258774-105259016:+:228|chr14:105258859-105258986:+:1193#chr14:105258855-105258982:-:916#chr14:105258810-105258944:-:274#chr14:105258806-105258944:+:323#chr14:105258778-105259019:-:47#chr14:105258774-105259016:+:145\n"
        variant_2 = "6455\t1\t1.5491866769945777E-4\tNC_000014.8\t105258893\tA\t1|6454|0|0\tG\tA\t0.9998450813323005\t0|0|1|0\t1915|1318|1320|1901\t0|0|0|0\t0|0|0|0\tchr14:105258855-105258982:-:1|chr14:105258859-105258986:+:2737#chr14:105258855-105258982:-:1929#chr14:105258810-105258944:-:595#chr14:105258806-105258944:+:686#chr14:105258778-105259019:-:95#chr14:105258774-105259016:+:374|0|0"
        variant_3 = "1706\t1601\t0.9384525205158265\tNC_000010.10\t123256344\tT\t104|1|0|1601\tT\tA\t0.9390029325513196\t0|52|52|0\t0|0|0|1\t0|0|0|0\t235|567|568|231\t0|chr10:123256234-123256344:+:1|0|chr10:123256234-123256344:+:416#chr10:123256234-123256346:-:588"
        variant_4 = "1148\t795\t0.6925087108013938\tNC_000010.10\t123256345\tT\t353|0|0|795\tT\tA\t0.6925087108013938\t0|182|171|0\t0|0|0|0\t0|0|0|0\t23|374|374|24\t0|0|0|chr10:123256234-123256344:+:1#chr10:123256234-123256346:-:587"
        variant_5 = "670\t567\t0.8462686567164179\tNC_000010.10\t123256346\tA\t567|0|0|103\tA\tT\t0.8462686567164179\t1|292|274|0\t0|0|0|0\t0|0|0|0\t0|57|45|1\tchr10:123256234-123256346:-:564|0|0|chr10:123256234-123256346:-:3"
        variant_6 = "1412\t700\t0.5042492917847026\tNC_000010.10\t123239112\tG\t712|700|0|0\tG\tA\t0.5042492917847026\t170|83|19|428\t213|72|25|402\t0|0|0|0\t0|0|0|0\tchr10:123239059-123239249:-:36#chr10:123239055-123239245:+:340#chr10:123238910-123239122:+:258#chr10:123238910-123239121:-:64|chr10:123239059-123239249:-:40#chr10:123239055-123239245:+:406#chr10:123238910-123239122:+:192#chr10:123238910-123239121:-:47|0|0"

        for line in [variant_1, variant_2, variant_3, variant_4, variant_5, variant_6]:
            alleles = extract_major_snv_allele('sample1', ref, line, self.nc_to_chr, 20, allele_filter)
            self.assertEqual(1, len(alleles.keys()))
            for key, info in alleles.items():
                (nc, pos, pos, ref_base, var_base) = key.split("#")
                try:
                    if int(pos) in variants[nc]:
                        variants[nc][int(pos)][key] = info
                    else:
                        variants[nc][int(pos)] = {key: info}
                except KeyError:
                    variants[nc] = {int(pos): {key: info}}
        self.maxDiff = None
        info = generate_multi_bp_variants("sample1", self.nc_to_chr,variants,20)
        self.assertEqual(4,len(info.keys()))
        self.assertEqual("14\t105258892\t105258893\tCA\tTG\tcomments: sample=sample1 variantAlleleRatio=0.4528856568157202 alleleFreq=3532,2927 readDepth=6463 Tumor_A=- Tumor_G=- Tumor_C=- Tumor_T=- Tumor_var_plusAmplicons=3 Tumor_var_minusAmplicons=3 Tumor_ref_plusAmplicons=3 Tumor_ref_minusAmplicons=3 Tumor_var_ampliconInfo=chr14:105258859-105258986:+:1193#chr14:105258855-105258982:-:916#chr14:105258810-105258944:-:274#chr14:105258806-105258944:+:323#chr14:105258778-105259019:-:47#chr14:105258774-105259016:+:145 Tumor_ref_ampliconInfo=chr14:105258859-105258986:+:1544#chr14:105258855-105258982:-:1010#chr14:105258810-105258944:-:324#chr14:105258806-105258944:+:363#chr14:105258778-105259019:-:52#chr14:105258774-105259016:+:228",info["NC_000014.8#105258892#105258893#CA#TG"])
        self.assertEqual("10\t123256344\t123256345\tTT\tAA\tcomments: sample=sample1 variantAlleleRatio=0.0609613130128957 alleleFreq=1601,104 readDepth=1706 Tumor_A=- Tumor_G=- Tumor_C=- Tumor_T=- Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=0 Tumor_ref_ampliconInfo=chr10:123256234-123256344:+:416#chr10:123256234-123256346:-:588",info["NC_000010.10#123256344#123256345#TT#AA"])
        self.assertEqual("10\t123256344\t123256346\tTTA\tAAT\tcomments: sample=sample1 variantAlleleRatio=0.0609613130128957 alleleFreq=1601,104 readDepth=1706 Tumor_A=- Tumor_G=- Tumor_C=- Tumor_T=- Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=1 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=0 Tumor_ref_ampliconInfo=chr10:123256234-123256344:+:416#chr10:123256234-123256346:-:588",info["NC_000010.10#123256344#123256346#TTA#AAT"])
        self.assertEqual("10\t123256345\t123256346\tTA\tAT\tcomments: sample=sample1 variantAlleleRatio=0.1537313432835821 alleleFreq=567,103 readDepth=670 Tumor_A=- Tumor_G=- Tumor_C=- Tumor_T=- Tumor_var_plusAmplicons=0 Tumor_var_minusAmplicons=0 Tumor_ref_plusAmplicons=0 Tumor_ref_minusAmplicons=1 Tumor_var_ampliconInfo=0 Tumor_ref_ampliconInfo=chr10:123256234-123256344:+:1#chr10:123256234-123256346:-:587",info["NC_000010.10#123256345#123256346#TA#AT"])

if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    unittest.main()
