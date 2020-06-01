#!/usr/bin/env python3

import pysam
from unittest import TestCase

import indeltools

reference = pysam.FastaFile("sample.fa")

class TestVariant(TestCase):
    def test_validate_allele(self):

        # None allele
        self.assertRaises(
            Exception, src.variant.Variant, "chrN", 12345, None, "C", reference
        )

        # empty allele
        self.assertRaises(
            Exception, src.variant.Variant, "chrN", 12345, "AA", "", reference
        )

        # non ACTG allele
        self.assertRaises(
            Exception, src.variant.Variant, "chrN", 12345, "T", "N", reference
        )

    def test_validate_locus(self):


        # No such a chromosome
        self.assertRaises(
            Exception, src.variant.Variant, "chrN", 12345, "A", "C", reference
        )

        # "N" in the genome
        #self.assertRaises(Exception, src.variant.Variant, "chr1", 1, "N", "N", REF_HG38)

        # " " (white-space) in the genome
        self.assertRaises(
            Exception, src.variant.Variant, "chrY", 999, " ", " ", REF_HG38
        )

    def setUp(self):
        self.var1 = src.variant.Variant("1", 14599, "T", "A", REF_HG37)

        self.var2 = src.variant.Variant("1", 14599, "T", "A", REF_HG37)
        self.var2.id = "."

        self.var3 = src.variant.Variant("1", 14599, "T", "A", REF_HG37)
        self.var3.id = "rs531646671"

        self.var4 = src.variant.Variant("1", 14599, "T", "A", REF_HG37)
        self.var4.id = "rs531646671"
        self.var4.id = "."
        self.var4.id = "rs707680"
        self.var4.is_common = True
        self.var4.is_common = False
        self.var4.is_common = True

        # self.var5 = src.variant.Variant("1", 14599, " ", " ", REF_HG37)

    def test_id(self):
        self.assertIsNone(self.var1.id)
        self.assertIsNone(self.var2.id)
        self.assertEqual(self.var3.id == "rs531646671", True)
        self.assertEqual(self.var4.id == "rs531646671,rs707680", True)
        self.assertEqual(self.var4.is_common == True, True)
        # self.assertEqual(self.var5.is_indel == False, True)


class TestSNV(TestCase):
    """Test SNV amd MNV cases
    """

    def setUp(self):
        # SNV
        self.snv1 = src.variant.Variant("chr1", 13148, "G", "A", REF_HG38)
        # non parsimonious expression
        self.snv2 = src.variant.Variant("chr1", 13146, "CCG", "CCA", REF_HG38)
        # MNV
        self.snv3 = src.variant.Variant("chr3", 45636145, "GCA", "TGC", REF_HG38)

        # non parsimonious expression
        self.snv4 = src.variant.Variant("chr3", 45636144, "TGCAA", "TTGCA", REF_HG38)
        # lower case input
        self.snv5 = src.variant.Variant("chr3", 45636145, "GcA", "tgC", REF_HG38)
        self.snv6 = src.variant.Variant("chr3", 45636144, "TgcAA", "ttgCa", REF_HG38)

    def test_equivalence(self):
        self.assertEqual(self.snv1 == self.snv2, True)
        self.assertEqual(self.snv3 == self.snv4, True)
        self.assertEqual(self.snv5 == self.snv6, True)


class TestIndel(TestCase):
    """Test indel cases
    """

    def setUp(self):

        # insertion
        self.ins1 = src.variant.Variant("chrX", 48791259, "G", "GAGCAC", REF_HG38)

        self.ins2 = src.variant.Variant("chrX", 48791265, "A", "AGCACA", REF_HG38)
        # non parsimonious expression
        self.ins3 = src.variant.Variant("chrX", 48791264, "CA", "CAGCACA", REF_HG38)
        # non parsimonious expression 2
        self.ins4 = src.variant.Variant("chrX", 48791264, "CAG", "CAGCACAG", REF_HG38)

        # deletion
        self.del1 = src.variant.Variant("4", 103827696, "CAAT", "C", REF_HG37)
        self.del2 = src.variant.Variant("4", 103827697, "AATA", "A", REF_HG37)
        # non parsimonious expression
        self.del3 = src.variant.Variant("4", 103827695, "ACAAT", "AC", REF_HG37)
        # non parsimonious expression 2
        self.del4 = src.variant.Variant("4", 103827695, "ACAATA", "ACA", REF_HG37)
        # lower case input
        self.ins5 = src.variant.Variant("chrX", 48791259, "g", "GAgcAc", REF_HG38)
        self.ins6 = src.variant.Variant("chrX", 48791265, "A", "aGCaca", REF_HG38)
        self.del5 = src.variant.Variant("4", 103827696, "cAAt", "c", REF_HG37)
        self.del6 = src.variant.Variant("4", 103827697, "AAtA", "A", REF_HG37)

    def test_equivalence(self):
        self.assertEqual(self.ins1 == self.ins2, True)
        self.assertEqual(self.ins1 == self.ins3, True)
        self.assertEqual(self.ins1 == self.ins4, True)

        self.assertEqual(self.del1 == self.del2, True)
        self.assertEqual(self.del1 == self.del3, True)
        self.assertEqual(self.del1 == self.del4, True)
        self.assertEqual(self.ins5 == self.ins6, True)
        self.assertEqual(self.del5 == self.del6, True)
