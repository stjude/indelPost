#!/usr/bin/env python3

import random
import numpy as np
#from .contig import Contig

from .pileup import (
    make_pileup,
    update_pileup,
    retarget,
    spliced_reference,
    update_read_info,
    check_overhangs,
    filter_spurious_overhangs,
)
from .equivalence import find_by_equivalence
from .softclip import find_by_softclip_split
from .localn import find_by_smith_waterman_realn, make_aligner
from .alleles import suggest_from_alignment
from .alleles_working import hard_phase_nearby_variants
from .utilities import get_local_reference, split


from indelpost.variant cimport Variant
from indelpost.contig cimport Contig 

from pysam.libcalignmentfile cimport AlignmentFile


random.seed(123)


cdef class VariantAlignment:
    cdef Variant target, __target
    cdef AlignmentFile bam
    cdef int window, retarget_window, mapqthresh
    cdef int downsamplethresh, basequalthresh, match_score
    cdef int mismatch_penalty, gap_open_penalty, gap_extension_penalty
    cdef float retarget_cutoff, __sample_factor
    cdef bint exclude_duplicates, retarget, perform_retarget, is_bad_overhang
    cdef list __pileup
    cdef Contig __contig
    
    def __cinit__(
        self,
        Variant target,
        AlignmentFile bam,
        int window=50,
        bint exclude_duplicates=True, 
        bint retarget=True,
        int retarget_window=30,
        float retarget_cutoff=0.6,
        int mapqthresh=20,
        int downsamplethresh=30000,
        int basequalthresh=20,
        int match_score=2,
        int mismatch_penalty=2,
        int gap_open_penalty=3,
        int gap_extension_penalty=1,
    ):
        
        if not target.is_non_complex_indel():
            raise Exception("Expected input is a non-complex indel")
        
        self.target = target
        self.__target = target.normalize()
        self.bam = bam
        self.window = window
        self.exclude_duplicates = exclude_duplicates
        self.perform_retarget = retarget
        self.retarget_window = retarget_window
        self.retarget_cutoff = retarget_cutoff
        self.mapqthresh = mapqthresh
        self.downsamplethresh = downsamplethresh  # should allow > 1000
        self.basequalthresh = basequalthresh
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_open_penalty = gap_open_penalty
        self.gap_extension_penalty = gap_extension_penalty

        self.__pileup, self.__contig = self.__parse_pileup()

    def __parse_pileup(self, contig=None, retargeted=False):
        """Dictize reads for target indel and make target indel template by consensus

        bam (pysam.AlignmentFile)
        target_indel (Variant): indel of interest 
        smith_waterman (bool): True to perform SW realn
        template (bool): None if no indel template supplied
        """

        # equivalence search
        if retargeted:
            pileup = self.__pileup
            sample_factor = self.__sample_factor
        else:
            pileup, self.__sample_factor = make_pileup(
                self.__target,
                self.bam,
                exclude_duplicates=self.exclude_duplicates,
                window=self.window,
                downsamplethresh=self.downsamplethresh,
                basequalthresh=self.basequalthresh,
            )

            self.__target, pileup = find_by_equivalence(
                self.__target,
                pileup,
                self.match_score,
                self.mismatch_penalty,
                self.gap_open_penalty,
                self.gap_extension_penalty,
            )

            contig = Contig(
                self.__target,
                preprocess_for_contig_construction(
                    self.__target,
                    pileup,
                    self.match_score,
                    self.mismatch_penalty,
                    self.gap_open_penalty,
                    self.gap_extension_penalty,
                ),
            )

            self.is_bad_overhang = False
            if contig.failed:

                within = min(self.retarget_window, len(self.__target.indel_seq) * 3)

                ans = check_overhangs(pileup)
                if ans:
                    intron, overhangs = ans[0], ans[1]
                    non_spurious_overhangs = filter_spurious_overhangs(
                        self.__target,
                        intron,
                        overhangs,
                        self.match_score,
                        self.mismatch_penalty,
                        self.gap_open_penalty,
                        self.gap_extension_penalty,
                    )
                    if not non_spurious_overhangs:
                        contig = Contig(self.__target, [])
                        self.is_bad_overhang = True
                        return pileup, contig
                    else:
                        res = retarget(
                            self.perform_retarget,
                            self.__target,
                            non_spurious_overhangs,
                            self.mapqthresh,
                            within,
                            self.retarget_cutoff,
                        )
                        if not res:
                            contig = Contig(self.__target, [])
                            self.is_bad_overhang = True
                            return pileup, contig
                else:
                    res = retarget(
                        self.perform_retarget,
                        self.__target,
                        pileup,
                        self.mapqthresh,
                        within,
                        self.retarget_cutoff,
                        self.match_score,
                        self.mismatch_penalty,
                        self.gap_open_penalty,
                        self.gap_extension_penalty,
                    )

                # if retargeted successfully -> make template based on the retarget
                if res:
                    self.__target, retarget_reads = res[0], res[1]
                    self.__target, self.__pileup = update_pileup(
                        pileup,
                        self.__target,
                        retarget_reads,
                        self.match_score,
                        self.mismatch_penalty,
                        self.gap_open_penalty,
                        self.gap_extension_penalty,
                    )

                    contig = Contig(
                        self.__target,
                        preprocess_for_contig_construction(
                            self.__target,
                            self.__pileup,
                            self.match_score,
                            self.mismatch_penalty,
                            self.gap_open_penalty,
                            self.gap_extension_penalty,
                        ),
                    )

                    # 2nd-pass using the retarget
                    return self.__parse_pileup(contig=contig, retargeted=True)

                # no target in this pileup
                else:
                    contig = Contig(self.__target, [])
                    return pileup, contig

        # soft-clip realn & SW realn
        if contig.qc_passed:
            pileup = find_by_softclip_split(self.__target, contig, pileup)
            pileup = find_by_smith_waterman_realn(
                self.__target,
                contig,
                pileup,
                self.match_score,
                self.mismatch_penalty,
                self.gap_open_penalty,
                self.gap_extension_penalty,
            )

        return pileup, contig

    def count_alleles(
        self, bint fwrv=False, bint by_fragment=False, int qualitywindow=0, int qualitythresh=0
    ):
        """Unique-count reads with and without target indel

        Args:
            pileup (list): list of dictized read (dict)
            target_indel (Variant)
        Returns:
        allele counts (tuple): (ref_count, alt_count)
        """

        cdef dict read
        
        pos, indel_len = self.target.pos, len(self.target.indel_seq)
        margin = min(indel_len / 2, 5) if indel_len > 3 else indel_len

        reads = self.__pileup
        if qualitywindow and qualitythresh:
            reads = [
                read
                for read in reads
                if is_quality_read(read, pos, qualitywindow, qualitythresh)
            ]

        fw_target = [
            read["read_name"]
            for read in reads
            if read["is_target"] and not read["is_reverse"]
        ]
        rv_target = [
            read["read_name"]
            for read in reads
            if read["is_target"] and read["is_reverse"]
        ]

        fw_non_target = [
            read["read_name"]
            for read in reads
            if count_as_non_target(read, pos, margin) and not read["is_reverse"]
        ]
        rv_non_target = [
            read["read_name"]
            for read in reads
            if count_as_non_target(read, pos, margin) and read["is_reverse"]
        ]

        fw_target = set(fw_target)
        rw_target = set(rv_target)
        fwrv_target = fw_target | rw_target

        fw_non_target = set(fw_non_target) - fwrv_target
        rv_non_target = set(rv_non_target) - fwrv_target
        if fwrv:
            return (
                (
                    int(len(fw_non_target) * self.__sample_factor),
                    int(len(rv_non_target) * self.__sample_factor),
                ),
                (
                    int(len(fw_target) * self.__sample_factor),
                    int(len(rv_target) * self.__sample_factor),
                ),
            )
        else:
            if by_fragment:
                fwrv_non_target = len(fw_non_target | rv_non_target)
                fwrv_target = len(fw_target | rw_target)
            else:
                fwrv_non_target = len(fw_non_target) + len(rv_non_target)
                fwrv_target = len(fw_target) + len(rw_target)

            return (
                int(fwrv_non_target * self.__sample_factor),
                int(fwrv_target * self.__sample_factor),
            )

    def suggest(self, dbsnp=None):
        return suggest_from_alignment(
            self.__target, self.__contig, self.__pileup, dbsnp
        )

    def to_complex(
        self,
        snv_neigborhood=15,
        indel_nearby=20,
        indel_repeat_thresh=2,
        mapq_thresh=20,
        dbsnp=None,
    ):
        return hard_phase_nearby_variants(
            self.__target,
            self.__contig,
            self.__pileup,
            snv_neigborhood,
            indel_nearby,
            indel_repeat_thresh,
            mapq_thresh,
            dbsnp,
        )


def is_quality_read(read, pos, qualitywindow, qualitythresh):

    try:
        lt_qual, rt_qual = read["lt_qual"], read["rt_qual"]
    except:
        # return True
        lt_qual, rt_qual = split(
            read["read_qual"],
            read["cigar_string"],
            pos,
            read["aln_start"],
            is_for_ref=False,
            reverse=False,
        )

    # return True
    if lt_qual and rt_qual:
        lt_median = np.median(lt_qual[-min(len(lt_qual), qualitywindow) :])
        rt_median = np.median(rt_qual[: min(len(rt_qual), qualitywindow)])

        return lt_median > qualitythresh and rt_median > qualitythresh


cdef bint count_as_non_target(dict read, int pos, int margin):
    if read["is_target"] or not read["is_covering"]:
        return False
    else:
        margin = 0 if read["is_spliced"] else margin
        covering_subread = read["covering_subread"]

        if covering_subread[0] + margin < pos < covering_subread[1] - margin:
            return True


cdef list preprocess_for_contig_construction(
    Variant target,
    list pileup,
    int match_score,
    int mismatch_penalty,
    int gap_open_penalty,
    int gap_extension_penalty,
):
    
    cdef dict read
    cdef int clips, nonclips

    if not pileup:
        return pileup

    targetpileup = [read for read in pileup if read["is_target"]]
    if not targetpileup:
        return targetpileup

    clipped_targetpileup = [
        read for read in targetpileup if "S" in read["cigar_string"]
    ]
    nonclipped_targetpileup = [
        read for read in targetpileup if not "S" in read["cigar_string"]
    ]

    clips = len(clipped_targetpileup)
    nonclips = len(nonclipped_targetpileup)

    if nonclips > 19:
        targetpileup = random.sample(nonclipped_targetpileup, 20)
    else:
        if nonclips + clips > 19:
            clipped_targetpileup = random.sample(clipped_targetpileup, 20 - nonclips)

        ref_seq, lt_len = get_local_reference(target, pileup)
        aligner = make_aligner(ref_seq, match_score, mismatch_penalty)
        ref_start = target.pos + 1 - lt_len
        
        is_gapped_aln = False
        clipped_targetpileup = [
            update_read_info(
                read,
                target,
                is_gapped_aln,
                gap_open_penalty,
                gap_extension_penalty,
                aligner,
                ref_seq,
                ref_start,
            )
            for read in clipped_targetpileup
            if not read.get("cigar_updated", False)
        ]

        targetpileup = nonclipped_targetpileup + clipped_targetpileup

    return targetpileup


#
#
# def suggest_del_seq(target_deletion, template, pileup):
#     chrom, pos, del_orig, genome = target_deletion.chrom, target_deletion.pos, target_deletion.indel_seq, target_deletion.genome
#
#
#     lt_ref, rt_ref = template["lt_ref"], template["rt_ref"]
#     lt_flank, rt_flank = template["lt_flank"], template["rt_flank"]
#
#     lt_shift, rt_shift = 0, 0
#     if has_mismatches(lt_ref):
#         lt_shift = pos_shift_by_exact_match(chrom, pos, lt_flank, genome, is_left=True)
#
#     if has_mismatches(rt_ref):
#         rt_shift = pos_shift_by_exact_match(chrom, pos, rt_flank, genome)
#
#     if lt_shift or rt_shift:
#         new_ref = genome.fetch(chrom, pos - 1 + lt_shift, pos + len(del_orig) + rt_shift)
#         new_alt = genome.fetch(chrom, pos - 1 + lt_shift, pos + lt_shift)
#         return pos + lt_shift, new_ref, new_alt
#     else:
#         return None
#
#
#
#
#
#
#
# def has_mismatches(flank_ref_seq):
#     res = False if flank_ref_seq.isupper() else True
#     return res
#
#
#
#
#
#
# def infer_indel_sequence(target_indel, pileup, second_pileup=None):
#     if skip_denovo_repeat_check(target_indel, pileup):
#         return target_indel.indel_seq
#
#     pos = target_indel.pos
#     del_seq = target_indel.indel_seq
#     denovo_repeats = [
#         check_for_denovo_repeat(d, pos, del_seq)
#         for d in pileup
#         if check_for_denovo_repeat(d, pos, del_seq)
#     ]
#     if denovo_repeats:
#         return most_common(denovo_repeats)
#     else:
#         return del_seq
#
#
# def skip_denovo_repeat_check(target_indel, pileup):
#     """Check if targe indel is a non-repetitive deletion
#        for polymorphism-induced de novo repeat checking
#     """
#     if target_indel.variant_type == "I":
#         return True
#
#     targets = [d for d in pileup if d["is_target"]]
#
#     if not targets:
#         return True
#
#     del_seq = target_indel.indel_seq
#     del_len = len(del_seq)
#     is_in_repeats = [
#         (del_seq == d["rt_flank"][:del_len]) or (del_seq == d["lt_flank"][-del_len:])
#         for d in targets
#     ]
#
#     if most_common(is_in_repeats):
#         return True
#
#     return False
#
#
# def check_for_denovo_repeat(read_dict, pos, del_seq):
#     if read_dict["is_target"] or not read_dict["is_covering"]:
#         return None
#
#     # do not infer from reads with lower MAPQ
#     if read_dict["mapq"] < 20:
#         return None
#
#     read_seq, read_pos, cigar_list = (
#         read_dict["read"].query_sequence,
#         read_dict["read_start"],
#         read_dict["cigar_list"],
#     )
#
#     read_idx = 0
#     prev_event = ""
#     for token in cigar_list:
#         event, event_len = token[-1], int(token[:-1])
#
#         if read_pos < pos:
#             read_pos = read_pos if event == "I" else (read_pos + event_len)
#         else:
#             break
#
#         prev_event = event
#         read_idx = read_idx if event == "D" else read_idx + event_len
#
#     # do not infer from softclipped bases
#     if prev_event == "S" or event == "S":
#         return None
#
#     # do not infer from read ends
#     if read_idx / len(read_seq) < 0.15 or read_idx / len(read_seq) > 0.85:
#         return None
#
#     diff = pos - read_pos
#
#     indel_len = len(del_seq)
#     quals = read_dict["read"].query_qualities
#     lt_flank, mid_seq, mid_quals, rt_flank = (
#         read_seq[: read_idx + diff],
#         read_seq[read_idx + diff : read_idx + diff + indel_len],
#         quals[read_idx + diff : read_idx + diff + indel_len],
#         read_seq[read_idx + diff + indel_len :],
#     )
#
#     # if repetitive repalace
#     if lt_flank and mid_seq and rt_flank:
#         if np.median(mid_quals) > 24 and (
#             mid_seq == lt_flank[-indel_len:] or mid_seq == rt_flank[:indel_len]
#         ):
#             return mid_seq
#     else:
#         return None
#
#
# def pos_shift_by_exact_match(chrom, pos, flank_seq, genome, is_left=False):
#     if is_left:
#         ref_seq = genome.fetch(chrom, pos - 50, pos)
#         ref_seq = ref_seq[::-1]
#         flank_seq = flank_seq[::-1].upper()
#     else:
#         ref_seq = genome.fetch(chrom, pos, pos - 50)
#
#     res = ref_seq.find(flank_seq)
#
#     if res == -1:
#         return 0
#
#     pos_shift = - res if is_left else res
#
#     return pos_shift
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
