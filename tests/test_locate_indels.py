#!/usr/bin/env python3

from unittest import TestCase

import indelpost


class TestData(TestCase):
    def setUp(self):
        self.cigar_0 = "3M678N48M8D31M19S"   
        self.start_pos_0 = 55174041
        
        self.cigar_1 = "25M8D28M8I15M6472N24M"
        self.start_pos_1 = 55174745
        
        self.cigar_2 = "26M11I3M6I1M2I51M"
        self.start_pos_2 = 111884609

    def test_locate_indels(self):
        
        # del with splicing
        self.assertEqual(indelpost.locate_indels(self.cigar_0, self.start_pos_0), ([], [(55174769,8)]))

        # ins and del
        self.assertEqual(indelpost.locate_indels(self.cigar_1, self.start_pos_1), ([(55174805,8)], [(55174769,8)]))

        # 3 ins
        self.assertEqual(indelpost.locate_indels(self.cigar_2, self.start_pos_2), ([(111884634, 11), (111884637, 6), (111884638, 2)], [])) 
