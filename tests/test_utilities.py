#!/usr/bin/env python3

from unittest import TestCase

import indelpost


class TestSplit(TestCase):
    def setUp(self):
        self.seq_0 = "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGAGAAGCAACATCTCCGAAAGCCAACAAGCCAACAAGGAAATCCTCGATGAAGC"
        self.seq_0_ref = "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAACAAGGAAATCCTCGATGAAGC"
        self.cigar_0 = "45M8D28M8I15M6472N5M"
        self.target_pos_00 = 55174805
        self.target_pos_01 = 55174769
        self.start_pos_0 = 55174725
        self.end_pos_0 = 55181297

        self.seq_1 = "GGTGCACACGGCTTGAGATGCCTGACGACTTTCGAGGAACGGAAGGCGGTTTACACCTTTGTGCTGAAGGTGAGTGACAAGGCTTTTCACACCCTGGGGC"
        self.seq_1_ref = "GGTGCACACGGCTTGAGATGCCTGACAACCTTTACACCTTTGTGCTGAAGGTGAGTGACAAGGCTTTTCACACCCTGGGGC"
        self.cigar_1 = "26M11I3M6I1M2I51M"
        self.target_pos_10 = 111884637
        self.start_pos_1 = 111884609
        self.end_pos_1 = 111884689

        self.seq_2 = "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAAATGAGGAATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT"
        self.seq_2_ref = "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAAGATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT"
        self.cigar_2 = "107M2I1M4I36M"
        self.target_pos_20 = 104850748
        self.target_pos_21 = 104850749
        self.start_pos_2 = 104850642
        self.end_pos_2 = 104850785

    def test_split(self):

        # forward for read
        self.assertEqual(
            indelpost.split(
                self.seq_0,
                self.cigar_0,
                self.target_pos_00,
                self.start_pos_0,
                is_for_ref=False,
                reverse=False,
            ),
            (
                "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGAGAAGCAACATCTCCGAAAGCCAAC",
                "AAGCCAACAAGGAAATCCTCGATGAAGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_0,
                self.cigar_0,
                self.target_pos_01,
                self.start_pos_0,
                is_for_ref=False,
                reverse=False,
            ),
            (
                "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATC",
                "AAGAGAAGCAACATCTCCGAAAGCCAACAAGCCAACAAGGAAATCCTCGATGAAGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_1,
                self.cigar_1,
                self.target_pos_10,
                self.start_pos_1,
                is_for_ref=False,
                reverse=False,
            ),
            (
                "GGTGCACACGGCTTGAGATGCCTGACGACTTTCGAGGAAC",
                "GGAAGGCGGTTTACACCTTTGTGCTGAAGGTGAGTGACAAGGCTTTTCACACCCTGGGGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_2,
                self.cigar_2,
                self.target_pos_20,
                self.start_pos_2,
                is_for_ref=False,
                reverse=False,
            ),
            (
                "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAA",
                "ATGAGGAATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_2,
                self.cigar_2,
                self.target_pos_21,
                self.start_pos_2,
                is_for_ref=False,
                reverse=False,
            ),
            (
                "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAAATG",
                "AGGAATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT",
            ),
        )

        # reverse for read
        self.assertEqual(
            indelpost.split(
                self.seq_0,
                self.cigar_0,
                self.target_pos_00,
                self.end_pos_0,
                is_for_ref=False,
                reverse=True,
            ),
            (
                "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGAGAAGCAACATCTCCGAAAGCCAAC",
                "AAGCCAACAAGGAAATCCTCGATGAAGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_0,
                self.cigar_0,
                self.target_pos_01,
                self.end_pos_0,
                is_for_ref=False,
                reverse=True,
            ),
            (
                "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATC",
                "AAGAGAAGCAACATCTCCGAAAGCCAACAAGCCAACAAGGAAATCCTCGATGAAGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_1,
                self.cigar_1,
                self.target_pos_10,
                self.end_pos_1,
                is_for_ref=False,
                reverse=True,
            ),
            (
                "GGTGCACACGGCTTGAGATGCCTGACGACTTTCGAGGAAC",
                "GGAAGGCGGTTTACACCTTTGTGCTGAAGGTGAGTGACAAGGCTTTTCACACCCTGGGGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_2,
                self.cigar_2,
                self.target_pos_20,
                self.end_pos_2,
                is_for_ref=False,
                reverse=True,
            ),
            (
                "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAA",
                "ATGAGGAATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_2,
                self.cigar_2,
                self.target_pos_21,
                self.end_pos_2,
                is_for_ref=False,
                reverse=True,
            ),
            (
                "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAAATG",
                "AGGAATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT",
            ),
        )

        # forward for ref
        self.assertEqual(
            indelpost.split(
                self.seq_0_ref,
                self.cigar_0,
                self.target_pos_00,
                self.start_pos_0,
                is_for_ref=True,
                reverse=False,
            ),
            (
                "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAAC",
                "AAGGAAATCCTCGATGAAGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_0_ref,
                self.cigar_0,
                self.target_pos_01,
                self.start_pos_0,
                is_for_ref=True,
                reverse=False,
            ),
            (
                "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATC",
                "AAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAACAAGGAAATCCTCGATGAAGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_1_ref,
                self.cigar_1,
                self.target_pos_10,
                self.start_pos_1,
                is_for_ref=True,
                reverse=False,
            ),
            (
                "GGTGCACACGGCTTGAGATGCCTGACAAC",
                "CTTTACACCTTTGTGCTGAAGGTGAGTGACAAGGCTTTTCACACCCTGGGGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_2_ref,
                self.cigar_2,
                self.target_pos_20,
                self.start_pos_2,
                is_for_ref=True,
                reverse=False,
            ),
            (
                "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAA",
                "GATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_2_ref,
                self.cigar_2,
                self.target_pos_21,
                self.start_pos_2,
                is_for_ref=True,
                reverse=False,
            ),
            (
                "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAAG",
                "ATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT",
            ),
        )

        # reverse for ref
        self.assertEqual(
            indelpost.split(
                self.seq_0_ref,
                self.cigar_0,
                self.target_pos_00,
                self.end_pos_0,
                is_for_ref=True,
                reverse=True,
            ),
            (
                "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAAC",
                "AAGGAAATCCTCGATGAAGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_0_ref,
                self.cigar_0,
                self.target_pos_01,
                self.end_pos_0,
                is_for_ref=True,
                reverse=True,
            ),
            (
                "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATC",
                "AAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAACAAGGAAATCCTCGATGAAGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_1_ref,
                self.cigar_1,
                self.target_pos_10,
                self.end_pos_1,
                is_for_ref=True,
                reverse=True,
            ),
            (
                "GGTGCACACGGCTTGAGATGCCTGACAAC",
                "CTTTACACCTTTGTGCTGAAGGTGAGTGACAAGGCTTTTCACACCCTGGGGC",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_2_ref,
                self.cigar_2,
                self.target_pos_20,
                self.end_pos_2,
                is_for_ref=True,
                reverse=True,
            ),
            (
                "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAA",
                "GATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT",
            ),
        )

        self.assertEqual(
            indelpost.split(
                self.seq_2_ref,
                self.cigar_2,
                self.target_pos_21,
                self.end_pos_2,
                is_for_ref=True,
                reverse=True,
            ),
            (
                "AAAGCATTTAAACTTATTTCCCAAGTTTTCTCCCAAATCCCATTTGGATACCTTAATACGTCTCTGGATGGAACTGATGTCTGGACGCTCATTGCTACTGCTGTCAAG",
                "ATGCCTAAAGATAAAAGAAAAACTCAGTGAAAATCT",
            ),
        )


class TestSpliced(TestCase):
    def setUp(self):
        self.cigar_0 = "3M678N48M8D31M19S"
        self.start_0 = 55174041
        self.end_0 = 55174827

        self.cigar_1 = "25M8D28M8I15M6472N24M"
        self.start_1 = 55174745
        self.end_1 = 55181316

    def test_spliced(self):
        self.assertEqual(
            indelpost.get_spliced_subreads(self.cigar_0, self.start_0, self.end_0),
            ([[55174041, 55174043], [55174722, 55174827]]),
        )
        self.assertEqual(
            indelpost.get_spliced_subreads(self.cigar_1, self.start_1, self.end_1),
            ([[55174745, 55174820], [55181293, 55181316]]),
        )


class TestCigar(TestCase):
    def setUp(self):
        self.cigar_0 = "1S17M678N48M8D31M4S"
        self.target_pos_0 = 1189
        self.target_pos_1 = 1221
        self.target_pos_2 = 1232
        self.target_pos_3 = 1234
        self.read_start_0 = 446

    def test_split_cigar(self):
        self.assertEqual(
            indelpost.split_cigar(self.cigar_0, self.target_pos_0, self.read_start_0),
            (["1S", "17M", "678N", "48M"], ["8D", "31M", "4S"]),
        )
        self.assertEqual(
            indelpost.split_cigar(self.cigar_0, self.target_pos_1, self.read_start_0),
            (["1S", "17M", "678N", "48M", "8D", "24M"], ["7M", "4S"]),
        )
        self.assertEqual(
            indelpost.split_cigar(self.cigar_0, self.target_pos_2, self.read_start_0),
            (["1S", "17M", "678N", "48M", "8D", "31M", "4S"], []),
        )
        self.assertEqual(
            indelpost.split_cigar(self.cigar_0, self.target_pos_3, self.read_start_0),
            None,
        )
