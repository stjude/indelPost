#!/usr/bin/env python3

import os
import pysam
from unittest import TestCase

import indelpost
print(indelpost.__file__)

tests = os.path.dirname(os.path.realpath(__file__))
reference = pysam.FastaFile(os.path.join(tests, "test.fa"))


class TestVariant(TestCase):
    def setUp(self):
        self.ins_0 = indelpost.Variant("chrZ", 1217, "A", "AAAGCCAAC", reference)
        self.ins_1 = indelpost.Variant("chrZ", 1225, "C", "CAAGCCAAC", reference)
        self.ins_2 = indelpost.Variant("chrZ", 1228, "G", "GCCAACAAG", reference)
        self.ins_3 = indelpost.Variant("chrZ", 1223, "A", "ACCAAGCCA", reference)

        self.del_0 = indelpost.Variant("chrZ", 454, "GGT", "G", reference)
        self.del_1 = indelpost.Variant("chrZ", 456, "TGT", "T", reference)
        self.del_2 = indelpost.Variant("chrZ", 453, "CGGT", "C", reference)
        self.del_3 = indelpost.Variant("chrZ", 454, "GGTG", "G", reference)

    def test_validate(self):

        # None allele
        self.assertRaises(
            Exception, indelpost.Variant, "chrZ", 400, "G", None, reference
        )

        # empty allele
        self.assertRaises(
            Exception, indelpost.Variant, "chrZ", 420, "AG", "", reference
        )

        # non ACTGN(actgn) allele
        self.assertRaises(
            Exception, indelpost.Variant, "chrZ", 460, "t", "B", reference
        )

        # No such a chromosome
        self.assertRaises(
            Exception, indelpost.Variant, "chrXYZ", 450, "G", "C", reference
        )

        # REF == ALT
        self.assertRaises(
            Exception, indelpost.Variant, "chrZ", 431, "C", "C", reference
        )

    def test_equivalence(self):
        self.assertEqual(self.ins_0 == self.ins_1 == self.ins_2, True)
        self.assertEqual(self.ins_0 == self.ins_3, False)
        self.assertEqual(self.del_0 == self.del_1, True)
        self.assertEqual(self.del_2 == self.del_3, True)

    def test_flank(self):
        self.assertEqual(self.ins_1.left_flank(10) == "GAAAGCCAAC", True)
        self.assertEqual(self.ins_1.right_flank(10) == "AAGGAAATCC", True)
        self.assertEqual(self.del_1.left_flank(10) == "TCGGCACGGT", True)
        self.assertEqual(self.del_1.right_flank(10) == "ATAAGGTAAG", True)

if __name__ == "__main__":
    unittest.main()
