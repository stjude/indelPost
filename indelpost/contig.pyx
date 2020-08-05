#cython: profile=True

import random
import numpy as np
from collections import OrderedDict

from .utilities import *

from indelpost.variant cimport Variant

from .consensus import make_consensus


random.seed(123)


cdef class Contig:
    def __cinit__(self, Variant target, list pileup, double low_consensus_thresh=0.7, int donwsample_lim=100):
        self.target = target
        self.pileup = pileup

        self.targetpileup = self.__preprocess(donwsample_lim)

        if self.targetpileup:
            self.__make_contig()
            self.failed = False
        else:
            self.qc_passed = False
            self.failed = True


    def __preprocess(self, donwsample_lim):
        targetpileup = [read for read in self.pileup if read["is_target"]]
        self.mapq = 0

        if not targetpileup:
            return targetpileup

        if len(targetpileup) > donwsample_lim:
            targetpileup = random.sample(targetpileup, donwsample_lim)

        self.mapq = np.median([read["mapq"] for read in targetpileup])

        return targetpileup


    def __make_contig(self):
        lt_consensus, rt_consensus = make_consensus(self.target, self.targetpileup)
        
        self.__index_by_genome_coord(lt_consensus[0], rt_consensus[0])

        self.start = lt_consensus[1]
        self.lt_reference_seq = lt_consensus[2]
        self.lt_consensus_seq = lt_consensus[3]
        self.lt_consensus_scores = lt_consensus[4]

        self.indel_seq = self.target.indel_seq

        self.rt_reference_seq = rt_consensus[2]
        self.rt_consensus_seq = rt_consensus[3]
        self.rt_consensus_scores = rt_consensus[4]
        self.end = rt_consensus[1]

        self.__profile_non_target_variants()

        self.qc_passed = self.__qc()


    def __index_by_genome_coord(self, lt_index, rt_index):
        self.lt_genomic_index = lt_index
        self.rt_genomic_index = rt_index

        tar = {self.target.pos: (self.target.ref, self.target.alt, 1.0)}

        genome_indexed_contig = rt_index
        genome_indexed_contig.update(lt_index)
        genome_indexed_contig.update(tar)
        self.contig_dict = OrderedDict(sorted(genome_indexed_contig.items()))


    def __profile_non_target_variants(self):
        non_target_variants = [
            Variant(self.target.chrom, k, v[0], v[1], self.target.reference)
            for k, v in self.contig_dict.items()
            if v[0] and v[0] != v[1] and k != self.target.pos
        ]
        self.non_target_indels = [var for var in non_target_variants if var.is_indel]
        self.mismatches = [var for var in non_target_variants if not var.is_indel]

        self.gaps = [
            str(len(var.indel_seq)) + var.variant_type for var in self.non_target_indels
        ]
        self.gaps.append(str(len(self.target.indel_seq)) + self.target.variant_type)

    def __qc(self):

        lt_n, lt_len = self.lt_consensus_seq.count("N"), len(self.lt_consensus_seq)
        rt_n, rt_len = self.rt_consensus_seq.count("N"), len(self.rt_consensus_seq)
        
        qc_stats ={}
        
        qc_stats["low_qual_base_frac"] = low_qual_fraction(self.targetpileup)
        
        qc_stats["clip_rate"] = sum(True for k, v in self.contig_dict.items() if not v[0]) / len(self.contig_dict)
        
        lt_n_proportion = lt_n / lt_len
        rt_n_proportion = rt_n / rt_len
        qc_stats["n_rate"] = (lt_n + rt_n) / (lt_len + rt_len)

        low_consensus_rate_lt = (
            sum(score < self.low_consensus_thresh for score in self.lt_consensus_scores) / lt_len
        )
        low_consensus_rate_rt = (
            sum(score < self.low_consensus_thresh for score in self.rt_consensus_scores) / rt_len
        )
       
        qc_stats["low_consensus_rate"] = (low_consensus_rate_lt * lt_len + low_consensus_rate_rt * rt_len) / (lt_len + rt_len)
        
        self.qc_stats = qc_stats

        if qc_stats["low_qual_base_frac"] > 0.2:
            return False
        elif qc_stats["clip_rate"] > 0.1:
            return False
        elif qc_stats["n_rate"] > 0.1:
            return False
        elif low_consensus_rate_lt > 0.2 or low_consensus_rate_rt > 0.2:
            return False
        else:
            return True

    
    def get_reference_seq(self, split=False):
        if self.failed:
            return None
        
        if split:
            if self.target.is_del:
                return self.lt_reference_seq, self.indel_seq, self.rt_reference_seq
            else:
                return  self.lt_reference_seq, "", self.rt_reference_seq
        else:
            refseq = (
                self.lt_reference_seq + self.indel_seq + self.rt_reference_seq
                if self.target.is_del
                else self.lt_reference_seq + self.rt_reference_seq
            )

            return refseq


    def get_contig_seq(self, split=False):
        if self.failed:
            return None
        
        if split:
            if self.target.is_ins:
                return self.lt_consensus_seq, self.indel_seq, self.rt_consensus_seq
            else:
                return self.lt_consensus_seq, "", self.rt_consensus_seq
        else:
            conseq = (
                self.lt_consensus_seq + self.indel_seq + self.rt_consensus_seq
                if self.target.is_ins
                else self.lt_consensus_seq + self.rt_consensus_seq
            )

            return conseq
