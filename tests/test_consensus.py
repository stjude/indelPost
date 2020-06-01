#!/usr/bin/env python3

import numpy as np
from unittest import TestCase

import indelpost


class TestData(TestCase):
    def setUp(self):
        self.refseq_lst_0 = ["ATTC", "ATTCG", "ATTCGG"]
        self.refseq_lst_1 = ["ATTN", "ATTNG", "ATTNGN"]
        self.refseq_lst_2 = ["ANCTTGC", "GTGC", "ATGA"]

    def test_consensus(self):
        res_0 = indelpost.consensus_refseq(self.refseq_lst_0)
        self.assertEqual(res_0[0], "ATTCGG")
        np.testing.assert_almost_equal(res_0[1], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

        res_1 = indelpost.consensus_refseq(self.refseq_lst_1)
        self.assertEqual(res_1[0], "ATTNGN")
        np.testing.assert_almost_equal(res_1[1], [1.0, 1.0, 1.0, 0.0, 1.0, 0.0])

        res_2 = indelpost.consensus_refseq(self.refseq_lst_2, left=True)
        self.assertEqual(res_2[0], "ANCATGC")
        np.testing.assert_almost_equal(
            res_2[1], [1.0, 0.0, 1.0, 0.3333333333333333, 1.0, 1.0, 0.6666666666666666]
        )

    def test_is_almost_same(self):
        self.assertEqual(
            indelpost.is_almost_same(
                "AATACAATGAGAC",
                "AATACAATGAGAC",
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            ),
            True,
        )
        self.assertEqual(
            indelpost.is_almost_same(
                "ACAGCAATGAGAC", "AATACAATGAGAC", [1, 0.5, 1, 0.9, 1, 1, 1, 1, 1, 1, 1]
            ),
            False,
        )
        self.assertEqual(
            indelpost.is_almost_same(
                "AATACAAAGAGAC", "AATTCAATGAGAC", [1, 1, 1, 0.5, 1, 1, 1, 1, 1, 1, 1]
            ),
            True,
        )
        self.assertEqual(indelpost.is_almost_same("AATA", "ACTA", [1, 1, 1, 1]), False)
