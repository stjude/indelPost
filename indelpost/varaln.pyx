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

from .alleles import phase_nearby_variants

from .utilities import get_local_reference, relative_aln_pos, split_cigar

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
    cdef readonly AlignmentFile bam
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
        int auto_adjust_extension_penalty_thresh=20,
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

        if (auto_adjust_extension_penalty_thresh <= 0 or 
            len(target.indel_seq) <= auto_adjust_extension_penalty_thresh
        ):    
            self.gap_extension_penalty = gap_extension_penalty
        else:
            self.gap_extension_penalty = 0

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

                within = self.retarget_window
                
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
                            self.window,
                            self.mapqthresh,
                            within,
                            self.retarget_cutoff,
                            self.match_score,
                            self.mismatch_penalty,
                            self.gap_open_penalty,
                            self.gap_extension_penalty,
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
    
    def __eq__(self, other):
        condition1 = (self.target_indel() == other.target_indel())
        condition2 = (self.bam.filename == other.bam.filename)
        return any((condition1, condition2))

    def __hash__(self):
        hashable = (self.target_indel(), self.bam.filename)
        return hash(hashable)

    def target_indel(self):
        return self.__target

    
    def fetch_reads(self, str how="target"):
        
        if how == "target":
            return [read["read"] for read in self.__pileup if read["is_target"]]
        elif how == "non_target":
            pos, indel_len = self.target.pos, len(self.target.indel_seq)
            margin = min(indel_len / 2, 5) if indel_len > 3 else indel_len
            return [read["read"] for read in self.__pileup if count_as_non_target(read, pos, margin)]
        elif how == "covering":
            return [read["read"] for read in self.__pileup if read["is_covering"]]
        else:
            return None       
            

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
        snv_neighborhood=15,
        indel_neighborhood=15,
        indel_repeat_thresh=5,
        sequence_complexity_thresh=0.01,
        dbsnp=None,
    ):
        return phase_nearby_variants(
            self.__target,
            self.contig,
            self.__pileup,
            self.mapqthresh,
            self.low_qual_frac_thresh,
            self.basequalthresh,
            snv_neighborhood,
            indel_neighborhood,
            indel_repeat_thresh,
            sequence_complexity_thresh,
            dbsnp,
        )
     
    def hard_phase(self):
        return phase_nearby_variants(
            self.__target,
            self.contig,
            self.__pileup,
            self.mapqthresh,
            self.low_qual_frac_thresh,
            self.basequalthresh,
            snv_neighborhood=15,
            indel_neighborhood=15,
            indel_repeat_thresh=5,
            sequence_complexity_thresh=0.02,
            dbsnp=None,
            hard=True,
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

    targetpileup = [read for read in pileup if read["is_target"] and not read["is_dirty"]]
    
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
        targetpileup = [right_aligner(read, target) for read in targetpileup]
    else:
        targetpileup = sorted(targetpileup, key=partial(centrality, target_pos=target.pos))
        
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
        
        targetpileup = [read for read in targetpileup if read is not None and (read.get("lt_cigar", None) and read.get("rt_cigar", None))]

        _targetpileup = [
            read for read in targetpileup if read.get("cigar_updated", False) 
        ]
        
        if len(_targetpileup) > 2:
            _targetpileup = _targetpileup[: min(20, int(len(_targetpileup) / 1.25 ))]
        
        if _targetpileup:    
            targetpileup = _targetpileup
        else:
            return targetpileup

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
    
    read = update_read_info(
        read, 
        target, 
        is_gapped_aln, 
        gap_open_penalty, 
        gap_extension_penalty, 
        aligner, 
        ref_seq, 
        ref_start
    )

    return right_aligner(read, target)


def right_aligner(read, target):
    """Right align indels around splice site"""
    
    if (
        "N" not in read["cigar_string"] 
        or (
            "I" not in read["cigar_string"] 
            and "D" not in read["cigar_string"]
        )
    ):
        return read     
   
    cigar_lst = read["cigar_list"]
    
    query_pos = 0
    ref_pos = read["read_start"]
    new_cigar = []
    prev_event = "A"
    skip_next = False
    for i, c in enumerate(cigar_lst):
        event, event_len = c[-1], int(c[:-1])
        
        if event_len < 0:
            return None

        query_move = 0 if event in ("D", "N", "H", "P") else event_len
        ref_move = 0 if event in ("I", "H", "P") else event_len

        if event in ("I", "D") and prev_event == "N":
            try:
                nxt_c = cigar_lst[i + 1]
                nxt_event, nxt_event_len = nxt_c[-1], int(nxt_c[:-1])
                if nxt_event != "M":
                    raise Exception
            except: 
                return None

            chrom, reference = target.chrom, target.reference
            padding_base = reference.fetch(chrom, ref_pos - 2, ref_pos - 1)
            if event == "I":
                ins_seq = read["read_seq"][query_pos: query_pos + event_len]
                ref = padding_base
                alt = padding_base + ins_seq
            else:
                del_seq = reference.fetch(chrom, ref_pos - 1, ref_pos - 1 + event_len)
                ref = padding_base + del_seq
                alt = padding_base
            
            right_aligned_vars = Variant(
                                    chrom, 
                                    ref_pos - 1, 
                                    ref, 
                                    alt, 
                                    reference
                                ).generate_equivalents()

            diff = max(v.pos for v in right_aligned_vars) - ref_pos + 1
            if diff > 0:
                new_cigar += [str(diff) + "M", str(event_len) + event, str(nxt_event_len - diff) + "M"]    
            else:
                return None

            ref_pos += query_move + nxt_event_len
            query_pos += ref_move + nxt_event_len
            skip_next = True

        else:
            if skip_next:
                skip_next = False
            else:
                query_pos += query_move
                ref_pos += ref_move
                new_cigar.append(c)
        
        prev_event = event
    
    read["cigar_list"] = new_cigar
    read["cigar_string"] = "".join(new_cigar)
    
    try:
        if target in right_aligned_vars:
            rt_aln_pos = target.pos + diff
            read["lt_cigar"], read["rt_cigar"] = split_cigar(read["cigar_string"], rt_aln_pos, read["read_start"])
            read["lt_flank"], read["rt_flank"] = split(
                                                      read["read_seq"], 
                                                      read["cigar_string"], 
                                                      rt_aln_pos, 
                                                      read["read_start"], 
                                                      is_for_ref=False, 
                                                      reverse=False
                                                 )
            read["lt_qual"], read["rt_qual"] = split(
                                                    read["read_qual"], 
                                                    read["cigar_string"], 
                                                    rt_aln_pos, 
                                                    read["read_start"], 
                                                    is_for_ref=False, 
                                                    reverse=False
                                               )
            read["lt_ref"], read["rt_ref"] = split(
                                                    read["ref_seq"], 
                                                    read["cigar_string"], 
                                                    rt_aln_pos, read["aln_start"], 
                                                    is_for_ref=True, 
                                                    reverse=False
                                             ) 
            read["target_right_shifted"] = rt_aln_pos
        else:
            read["lt_cigar"], read["rt_cigar"] = split_cigar(read["cigar_string"], target.pos, read["read_start"])
    except:
        pass

    return read
