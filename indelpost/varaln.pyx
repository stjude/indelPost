#cython: profile=False

cimport cython
import random
import numpy as np
from functools import partial

from .pileup import (
    update_pileup,
    retarget,
    update_read_info,
    check_overhangs,
    filter_spurious_overhangs,
)
from .gappedaln import find_by_normalization
from .softclip import find_by_softclip_split
from .localn import find_by_smith_waterman_realn, make_aligner

from .alleles import hard_phase_nearby_variants

from .utilities import get_local_reference, relative_aln_pos

from indelpost.pileup cimport make_pileup
from indelpost.utilities cimport split
from indelpost.contig cimport Contig 

from indelpost.variant cimport Variant

from pysam.libcalignmentfile cimport AlignmentFile

random.seed(123)

cdef class VariantAlignment:
    """Class to work with indels in the BAM file.

    Parameters
    ----------
    variant : Variant
        :class:`~indelpost.Variant` object representing the target indel.
    
    bam : pysam.AlignmentFile
        BAM file supplied as
        `pysam.AlignmentFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile>`__ object.

    window : integer
        indels will ba analyzed within the input indel position +/- window. Default to 50.

    exclude_duplicates : bool
        True to exclude reads with the 
        `MarkDuplicate <https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates>`__ signature.
        Default to True. 

    retarget : bool
        True to re-target the target if the input indel is not found after left-aligning the pileup.
        Default to True
    
    retarget_window : integer
        retarggeting will be performed within the input indel position +/- retarget_window. Default to 30.
    """
    cdef Variant target, __target
    cdef AlignmentFile bam
    cdef int window, retarget_window, mapqthresh
    cdef int downsamplethresh, basequalthresh, match_score
    cdef int mismatch_penalty, gap_open_penalty, gap_extension_penalty
    cdef float retarget_cutoff, __sample_factor, low_qual_frac_thresh
    cdef bint exclude_duplicates, retarget, perform_retarget
    cdef list __pileup
    cdef readonly is_spurious_overhang
    cdef readonly Contig contig
    
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
        float low_qual_frac_thresh=0.2,
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
        self.low_qual_frac_thresh = low_qual_frac_thresh
        self.downsamplethresh = downsamplethresh  # should allow > 1000
        self.basequalthresh = basequalthresh
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_open_penalty = gap_open_penalty
        self.gap_extension_penalty = gap_extension_penalty

        self.__pileup, self.contig = self.__parse_pileup()
    
    cdef __parse_pileup(self, Contig contig=None, bint retargeted=False):
        """Dictize reads for target indel and make target indel template by consensus

        bam (pysam.AlignmentFile)
        target_indel (Variant): indel of interest 
        smith_waterman (bool): True to perform SW realn
        template (bool): None if no indel template supplied
        """

        cdef list pileup
        
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
            
            for read in pileup:
                if "H" in read["read"].cigarstring:
                    print(read)
            
            
            self.__target, pileup, exptension_penalty_used = find_by_normalization(
                self.__target,
                pileup,
                self.window,
                self.match_score,
                self.mismatch_penalty,
                self.gap_open_penalty,
                self.gap_extension_penalty,
                self.basequalthresh,
            )

            if self.target != self.__target:
                    self.__target, pileup = update_pileup(
                        pileup,
                        self.__target,
                        self.window,
                        self.match_score,
                        self.mismatch_penalty,
                        self.gap_open_penalty,
                        self.gap_extension_penalty,
                        self.basequalthresh,
                        bypass_search=True
                    )

            contig = Contig(
                self.__target,
                preprocess_for_contig_construction(
                    self.__target,
                    self.target,
                    pileup,
                    self.window,
                    self.match_score,
                    self.mismatch_penalty,
                    self.gap_open_penalty,
                    exptension_penalty_used,
                ),
                self.basequalthresh,
                self.mapqthresh,
            )

            self.is_spurious_overhang = False
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
                        contig = Contig(self.__target, [], self.basequalthresh, self.mapqthresh)
                        self.is_spurious_overhang = True
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
                            contig = Contig(self.__target, [], self.basequalthresh, self.mapqthresh)
                            self.is_spurious_overhang = True
                            return pileup, contig
                else:
                    res = retarget(
                        self.perform_retarget,
                        self.__target,
                        pileup,
                        self.window,
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
                        self.window,
                        self.match_score,
                        self.mismatch_penalty,
                        self.gap_open_penalty,
                        self.gap_extension_penalty,
                        self.basequalthresh,
                        bypass_search=True,
                    )
                    
                    contig = Contig(
                        self.__target,
                        preprocess_for_contig_construction(
                            self.__target,
                            self.__target,
                            self.__pileup,
                            self.window,
                            self.match_score,
                            self.mismatch_penalty,
                            self.gap_open_penalty,
                            exptension_penalty_used,
                        ),
                        self.basequalthresh,
                        self.mapqthresh
                    )

                    # 2nd-pass using the retarget
                    return self.__parse_pileup(contig=contig, retargeted=True)

                # no target in this pileup
                else:
                    
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
        
        
            contig = Contig(
                self.__target,
                preprocess_for_contig_construction(
                    self.__target,
                    self.target,
                    pileup,
                    self.window,
                    self.match_score,
                    self.mismatch_penalty,
                    self.gap_open_penalty,
                    self.gap_extension_penalty,
                ),
                self.basequalthresh,
                self.mapqthresh,
            )
        
        
        
        return pileup, contig


    def target_indel(self):
        return self.__target

    
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


    def to_complex(
        self,
        snv_neigborhood=15,
        indel_nearby=15,
        indel_repeat_thresh=5,
        sequence_complexity_thresh=0.02,
        dbsnp=None,
    ):
        return hard_phase_nearby_variants(
            self.__target,
            self.contig,
            self.__pileup,
            self.mapqthresh,
            self.low_qual_frac_thresh,
            self.basequalthresh,
            snv_neigborhood,
            indel_nearby,
            indel_repeat_thresh,
            sequence_complexity_thresh,
            dbsnp,
        )


def is_quality_read(read, pos, qualitywindow, qualitythresh):

    try:
        lt_qual, rt_qual = read["lt_qual"], read["rt_qual"]
    except:
        lt_qual, rt_qual = split(
            read["read_qual"],
            read["cigar_string"],
            pos,
            read["aln_start"],
            is_for_ref=False,
            reverse=False,
        )

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


def centrality(read, target_pos):
    relative_pos = relative_aln_pos(read["ref_seq"], read["cigar_list"], read["aln_start"], target_pos)
    return abs(0.5 - relative_pos)


cdef list preprocess_for_contig_construction(
    Variant target,
    Variant orig_target,
    list pileup,
    int window,
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
        and (read.get("lt_cigar", None) and read.get("rt_cigar", None)) 
    ]

    
    clips = len(clipped_targetpileup)
    nonclips = len(nonclipped_targetpileup)

    if target == orig_target and nonclips > 19:
        targetpileup = random.sample(nonclipped_targetpileup, 20)
    else:
        targetpileup = sorted(targetpileup, key=partial(centrality, target_pos=target.pos))
        
        if len(targetpileup) > 2:
            targetpileup = targetpileup[: min(20, int(len(targetpileup) / 1.25 ))]
        
        unspl_ref_seq, unspl_lt_len = get_local_reference(orig_target, pileup, window, unspliced=True) 
        unspl_aligner = make_aligner(unspl_ref_seq, match_score, mismatch_penalty)
        unspl_start = orig_target.pos + 1 - unspl_lt_len 

        is_gapped_aln = False
        targetpileup = [
            update_spliced_read_info(
                read,
                target,
                orig_target,
                is_gapped_aln,
                window,
                match_score,
                mismatch_penalty,
                gap_open_penalty,
                gap_extension_penalty,
            )
            if "N" in read["cigar_string"] else 
            update_read_info(
                read,
                target,
                is_gapped_aln,
                gap_open_penalty,
                gap_extension_penalty,
                unspl_aligner,
                unspl_ref_seq,
                unspl_start,
            ) 
            for read in targetpileup
        ]
        
        updated_cigars = [read for read in targetpileup if read.get("cigar_updated", False)]
        
        #for read in updated_cigars:
        #    if not read.get("lt_cigar", False):
        #        print(read)
        
        if not updated_cigars:    
            return targetpileup
        else:
            targetpileup = updated_cigars

    return targetpileup


def update_spliced_read_info(
    read, 
    target, 
    orig_target, 
    is_gapped_aln, 
    window, 
    match_score, 
    mismatch_penalty, 
    gap_open_penalty, 
    gap_extension_penalty
):
    ref_seq, lt_len = get_local_reference(orig_target, [read], window)
    aligner = make_aligner(ref_seq, match_score, mismatch_penalty)
    ref_start = orig_target.pos + 1 - lt_len
    return update_read_info(
        read, 
        target, 
        is_gapped_aln, 
        gap_open_penalty, 
        gap_extension_penalty, 
        aligner, 
        ref_seq, 
        ref_start
    )

