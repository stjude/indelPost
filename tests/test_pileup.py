#!/usr/bin/env python3

import os
import pysam
from unittest import TestCase

import indelpost

tests = os.path.dirname(os.path.realpath(__file__))

bam = pysam.AlignmentFile(os.path.join(tests, "test.bam"))
reference = pysam.FastaFile(os.path.join(tests, "test.fa"))


class TestData(TestCase):
    def setUp(self):
        self.cigar_0 = "3M678N48M8D31M19S"
        self.start_pos_0 = 55174041

        self.cigar_1 = "25M8D28M8I15M6472N24M"
        self.start_pos_1 = 55174745

        self.cigar_2 = "26M11I3M6I1M2I51M"
        self.start_pos_2 = 111884609

        self.var_3 = indelpost.Variant("chrZ", 1228, "G", "GCCAACAAG", reference)
        self.cigar_3 = "1S16M678N48M8D31M8I3M"
        self.start_pos_3 = 447

        self.var_4 = indelpost.Variant("chrZ", 1189, "C", "CAAGGAATT", reference)
        self.cigar_4 = "24M678N48M8D28M"
        self.start_pos_4 = 440

        self.var_5 = indelpost.Variant("chrZ", 1164, "AGA", "A", reference)
        self.cigar_5 = "24M678N23M2D23M8D28M"
        self.start_pos_5 = 440

    def test_locate_indels(self):

        # del with splicing
        self.assertEqual(
            indelpost.locate_indels(self.cigar_0, self.start_pos_0),
            ([], [(55174769, 8)]),
        )

        # ins and del
        self.assertEqual(
            indelpost.locate_indels(self.cigar_1, self.start_pos_1),
            ([(55174805, 8)], [(55174769, 8)]),
        )

        # 3 ins
        self.assertEqual(
            indelpost.locate_indels(self.cigar_2, self.start_pos_2),
            ([(111884634, 11), (111884637, 6), (111884638, 2)], []),
        )

    def test_leftalign_cigar(self):
        self.assertEqual(
            indelpost.leftalign_cigar(self.cigar_3, self.var_3, self.start_pos_3),
            "1S16M678N48M8D20M8I14M",
        )

        self.assertEqual(
            indelpost.leftalign_cigar(self.cigar_4, self.var_4, self.start_pos_4),
            self.cigar_4,
        )

        self.assertEqual(
            indelpost.leftalign_cigar(self.cigar_5, self.var_5, self.start_pos_5),
            "24M678N21M2D25M8D28M",
        )
