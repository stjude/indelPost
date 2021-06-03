#cython: embedsignature=True
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
from .contig import compare_contigs
from .utilities import get_local_reference, relative_aln_pos, split_cigar

from indelpost.pileup cimport make_pileup
from indelpost.utilities cimport split
from indelpost.contig cimport Contig, FailedContig 

from indelpost.variant cimport Variant, NullVariant

from pysam.libcalignmentfile cimport AlignmentFile

random.seed(123)

cdef class VariantAlignment:
    """This class accepts the target indel as :class:`~indelpost.Variant` 
    and the BAM file as `pysam.AlignmentFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile>`__ .
    Upon creation, the target indel is searched through realignment.
    Evaluation of the equality between :class:`~indelpost.VariantAlignment` objects triggers local-phasing around the indels 
    to test the alignment identity. 


    Parameters
    ----------
    variant : Variant
        :class:`~indelpost.Variant` object representing the target indel.
    
    bam : pysam.AlignmentFile
        BAM file supplied as
        `pysam.AlignmentFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile>`__ object.

    window : integer
        analyzes region input_indel_pos +/- window (defalt 50).

    exclude_duplicates : bool
        True (default) to exclude reads marked duplicate by
        `MarkDuplicate <https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates>`__. 
    
    downsample_threshold : integer
        downsamples to the threshold if the covearge at the input locus exceeds the threshold (default to 1000).
        
    mapping_quality_threshold : integer
        reads with a mapping quality (MAPQ) < the threshold (default 1) will be filtered. 
        
    base_quality_threshold : integer
        non-reference base-calls with a quality score < the threshold (default 20) will be masked as "N" in the analysis.

    retarget_search_window : integer
        attempts to re-target in input_indel_pos +/- window (defalt 30) if the exact match is not found after realignment.

    retarget_similarity_cutoff : float
        re-targets to the indel sequence with a highest 
        `Ratcliff/Obershelp similarity score <https://docs.python.org/3/library/difflib.html#difflib.get_close_matches>`__ > cutoff (default 0.7).
   
    match_score : integer
        score (default 3) for matched bases in Smith-Waterman local alignment.

    mismatch_penaly : interger
        penalty (default 2) for mismatched bases in Smith-Waterman local alignment.
    
    gap_open_penalty : integer
        penalty (default 3) to create a gap in Smith-Waterman local alignment.
    
    gap_extension_penalty : integer 
        penalty (default 1) to expend a created gap by one base in Smith-Waterman local alignment.

    auto_adjust_extension_penalty : bool
        True (default) to auto-adjust the gap open and extend penalties to find the input indel by realignment.
    """
    def __cinit__(
        self,
        Variant target,
        AlignmentFile bam,
        int window=50,
        bint exclude_duplicates=True, 
        int retarget_search_window=30,
        float retarget_similarity_cutoff=0.7,
        int mapping_quality_threshold=1,
        int downsample_threshold=1000,
        int base_quality_threshold=20,
        int match_score=3,
        int mismatch_penalty=2,
        int gap_open_penalty=3,
        int gap_extension_penalty=1,
        bint auto_adjust_extension_penalty=True,
    ):
        
        self.target = target
        
        if not target.is_non_complex_indel() and target.is_indel:
            
            if auto_adjust_extension_penalty:
                decomposed_variants = target.decompose_complex_variant()
            else:
                decomposed_variants = target.decompose_complex_variant(match_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty)
            
            decomposed_indels = [i for i in decomposed_variants if i.is_indel]
            self.__target = max(decomposed_indels, key=lambda l : len(l.indel_seq))
            self.target = self.__target
        else:
            self.__target = target.normalize()

        self.bam = bam
        self.window = window
        self.exclude_duplicates = exclude_duplicates
        self.retarget_window = retarget_search_window
        self.retarget_cutoff = retarget_similarity_cutoff
        self.mapqthresh = mapping_quality_threshold
        self.downsamplethresh = downsample_threshold 
        self.basequalthresh = base_quality_threshold
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_open_penalty = gap_open_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.auto_adjust_extension_penalty = auto_adjust_extension_penalty
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
            
            self.__target, pileup, exptension_penalty_used, self._observed_pos = find_by_normalization(
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
                
                grid = generate_grid(self.auto_adjust_extension_penalty, 
                                     self.gap_open_penalty,
                                     self.gap_extension_penalty,
                                     self.__target,
                       )
                
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
                        res = grid_search(
                            self.__target, 
                            non_spurious_overhangs,
                            self.window,
                            self.mapqthresh,
                            within,
                            self.retarget_cutoff,
                            self.match_score,
                            self.mismatch_penalty,
                            grid
                        )
                         
                        if res:
                            self.gap_open_penalty, self.gap_extension_penalty = res[2], res[3]
                        else:      
                            contig = Contig(self.__target, [], self.basequalthresh, self.mapqthresh)
                            self.is_spurious_overhang = True
                            return pileup, contig
                else:
                    res = grid_search(
                        self.__target,
                        pileup,
                        self.window,
                        self.mapqthresh,
                        within,
                        self.retarget_cutoff,
                        self.match_score,
                        self.mismatch_penalty,
                        grid,
                    )
                    
                    if res:
                        self.gap_open_penalty, self.gap_extension_penalty = res[2], res[3]
                        
                # if retargeted successfully -> make contig based on the retarget
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
                            self.gap_extension_penalty,
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
            
            orig_contig = contig
            
            # realign reads that are not through retarget path
            if not retargeted:
               cutoff = 1.0
               within = 30 
               
               target = [read for read in pileup if read["is_target"]]
               nontarget = [read for read in pileup if not read["is_target"]]
              
               grid = generate_grid(self.auto_adjust_extension_penalty, 
                                    self.gap_open_penalty,
                                    self.gap_extension_penalty,
                                    self.__target,
                      )
               
               res = grid_search(
                       self.__target,
                       nontarget,
                       self.window,
                       self.mapqthresh,
                       within,
                       cutoff,
                       self.match_score,
                       self.mismatch_penalty,
                       grid,
                   )
              
               if res:
                   nontarget = [read for read in nontarget if read not in res[1]]
                   
                   pileup = target + res[1] + nontarget
                   self.gap_open_penalty, self.gap_extension_penalty = res[2], res[3]
                    
                   self.__target, pileup = update_pileup(
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
                          
                   #equivalent but different position
                   if self.__target == res[0]:
                       self.__target = res[0]

               else:
                  pileup = target + nontarget    
                    
            pileup = find_by_softclip_split(self.__target, contig, pileup)
            
            pileup = find_by_smith_waterman_realn(
                self.__target,
                contig,
                pileup,
                self.match_score,
                self.mismatch_penalty,
                self.gap_open_penalty,
                self.gap_extension_penalty,
                self.basequalthresh
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

            contig = compare_contigs(orig_contig, contig, self.__target.pos)
        
        return pileup, contig
    

    def __eq__(self, other):
        # alignment equivalence
        my_contig, other_contig = self.contig, other.contig

        if my_contig.failed or other_contig.failed:
            return False

        # check eq in phased form
        phasing_mode = "local"
        my_phased, other_phased = self.phase(how=phasing_mode), other.phase(how=phasing_mode) 
         
        return (my_phased == other_phased)


    def __hash__(self):
        hashable = self.phase(how="local")
        return hash(hashable)
   

    def get_contig(self):
        """returns :class:`~indelpost.Contig` object built from reads supporting the target indel.
        :class:`~indelpost.FailedContig` is returned if contig assembly is not successful.
        """
        contig = self.contig
        if contig and not contig.failed:
            return contig
        else:
            failed = FailedContig()
            alt_cnt = self.count_alleles()[1]
            if alt_cnt:
          
                dirty_target_pileup = [read["is_dirty"] for read in self.__pileup if read["is_target"]]
                if sum(dirty_target_pileup) == len(dirty_target_pileup):
                    failed.is_low_quality = True  
                else:
                    failed.failed_anyway = True
            else:
                failed.target_not_found = True

            return failed
         
    
    def get_target_indel(self):
        """returns :class:`~indelpost.Variant` object representing the actual target indel analyzed. 
        The target indel may be different from the input indel due to the internal realignment. 
        :class:`~indelpost.NullVariant` is returned when no target is found.
        """
        alt_cnt = self.count_alleles()[1]
        if alt_cnt:
            return self.__target
        else:
            return NullVariant(self.__target.chrom, self.__target.pos, self.__target.reference)

    
    def fetch_reads(self, how='target'):
        """returns a `list <https://docs.python.org/3/library/stdtypes.html#list>`__ 
        of `pysam.AlignedSegment <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment>`__ objects, 
        depending on the strategy specified. 

        Parameters
        ----------
            how : string
                specifies the read fetch strategy. "target" (default) to retrieve reads supporting the target.
                "non_target" for reads w/o target. "covering" for reads covering the locus w/ and w/o the target.
        """
        
        if how == "target":
            return [read["read"] for read in self.__pileup if read["is_target"]]
        elif how == "non_target":
            pos, indel_len = self.target.pos, len(self.target.indel_seq)
            margin = min(indel_len / 2, 5) if indel_len > 3 else indel_len
            return [read["read"] for read in self.__pileup if count_as_non_target(read, pos, margin)]
        elif how == "covering":
            return [read["read"] for read in self.__pileup if read["is_covering"]]
        else:
            raise Exception("fetch stragety must be either of target, non_target, covering")
            

    def count_alleles(
        self, fwrv=False, by_fragment=False, estimated_count=False, quality_window=None, quality_threshold=None
    ):
        """returns a `tuple <https://docs.python.org/3/library/stdtypes.html#tuple>`__ of
        read counts:  (#ref reads, #alt reads).

        Parameters
        ----------
            fwrv : bool
                breaks down into the forward and reverse read counts.
                ( (ref-fw, ref-rv), (alt-fw, alt-rv) ).
            by_fragment : bool
                counts by fragment. Overlapping fw and rv reads are counted as one.
            estimated_count : bool
                True to return estimated count when the coverage is higher than :attr:`~indelpost.VariantAlignment.downsample_threshold`.
            quality_window : integer
                specifies the range of base call quality filter. indel pos +/- quality_window. 
            quality_threshold : integer
                filters reads with the median base call quality in indel pos +/- quality_window < quality_threshold.
        """

        cdef dict read
        
        pos, indel_len = self._observed_pos, len(self.target.indel_seq)
        margin = min(indel_len / 2, 5) if indel_len > 3 else indel_len

        reads = self.__pileup
        if quality_window and quality_threshold:
            reads = [
                read
                for read in reads
                if is_quality_read(read, pos, quality_window, quality_threshold)
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

        est = self.__sample_factor if estimated_count else 1

        if fwrv:
            return (
                (
                    int(len(fw_non_target) * est),
                    int(len(rv_non_target) * est),
                ),
                (
                    int(len(fw_target) * est),
                    int(len(rv_target) * est),
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
                int(fwrv_non_target * est),
                int(fwrv_target * est),
            )


    def phase(
        self,
        how="local",
        local_threshold=20,
        longest_common_substring_threshold=15,
        indel_repeat_threshold=10,
        mutation_density_threshold=0.05,
    ):
        """returns a :class:`~indelpost.Variant` object represeting a phased target indel. 
        :class:`~indelpost.NullVariant` is returned if no target is found.
        Phasing is limited to the exon where the target indel is located for RNA-Seq. 
        Refer to :meth:`~indelpost.Contig.get_phasables()` to retrieve phasable variants across exons.
        
        Parameters
        ----------
        how : string
            - "local" (default) phasing by recursive elimination of longest common substrings (LCS) and locality  scoring (see local_threshold).
            - "greedy" phase all phasable events.
            - "complex" phasing by the "local" option + exclusivity check. The exclusivity check removes phasable events that are also observed in non-target reads such as nearby homozygous germline events.
        local_threshold : integer
            local (and complex) phasing method checks if phasable SNVs are locally clustered around the indel. 
            For i-th base apart from indel, the score gains one if mismatch and, if match, loses 
            1*min(j/local_thresh, 1) where j = i for short indels (< 10-nt) and 0.6*i for longer indels.      
        longest_common_substring_threshold : integer
            removes common substrings between the reference and contig assembled from target reads that are longer than longest_common_substring_threshold (default 15). 
        indel_repeat_threshold : integer
            do not phase indels that repeat more than indel_repeat_threshold (default 10)
        mutation_density_threshold : float
            do not phase if the pileup contains too many mutations (possibly error-prone dirty region). 
            In non-target reads, #non-ref bases/#all bases > mutation_density_threshold (default 0.05)
        """
        if how == "complex":
            hard, to_complex = False, True
        elif how == "greedy":
            hard, to_complex = True, False
        elif how == "local":
            hard, to_complex = False, False
        else:
            raise Exception("phasing stragety must be either of local, greedy, complex")
        
        return phase_nearby_variants(
            self.__target,
            self.contig,
            self.__pileup,
            self.basequalthresh,
            local_threshold,
            longest_common_substring_threshold,
            indel_repeat_threshold,
            mutation_density_threshold,
            hard,
            to_complex
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

        if covering_subread[0] + margin <= pos <= covering_subread[1] - margin:
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
    
    if target == orig_target and nonclips > 9:
        random.seed(123)
        targetpileup = random.sample(nonclipped_targetpileup, 10)
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

            indel_len = len(target.indel_seq)
            if target.is_ins:
                read["rt_flank"] = read["rt_flank"][indel_len :]
                read["rt_qual"] = read["rt_qual"][indel_len :]
            else:
                read["rt_ref"] = read["rt_ref"][indel_len :]
        else:
            read["lt_cigar"], read["rt_cigar"] = split_cigar(read["cigar_string"], target.pos, read["read_start"])
    except:
        pass

    return read


def generate_grid (auto_adjust_extension_penalty,
                   gap_open_penalty,
                   gap_extension_penalty,
                   target,
):
    if auto_adjust_extension_penalty:
        if (gap_open_penalty, gap_extension_penalty) != (3, 1):
            if len(target.indel_seq) < 20:
                grid = [(gap_open_penalty, gap_extension_penalty),  
                            (3, 0), (4, 0), (5, 0), (3, 1), (4, 1), (5, 1)
                       ]
            else:
                grid = [(gap_open_penalty, gap_extension_penalty),
                            (3, 1), (3, 0), (5, 1), (5, 0), (4, 1), (5, 1)
                       ]
        else:
            if len(target.indel_seq) < 20:
                grid = [(3, 1), (3, 0), (5, 1), (5, 0), (4, 1), (4, 0)]
            else:
                grid = [(3, 0), (3, 1), (5, 1), (5, 0), (4, 1), (4, 0)] 
    else:
        grid = [(gap_open_penalty, gap_extension_penalty)]
    
    return grid


def grid_search(
    target, 
    pileup,
    window,
    mapq_thresh,
    within,
    retarget_cutoff,
    match_score,
    mismatch_penalty,
    grid
):
    # grid = [(gap.open, gap.ext)]
    h = 0
    responses, scores, hs = [], [], []
    while h < len(grid):
        res = retarget(
            target,
            pileup,
            window,
            mapq_thresh,
            within,
            retarget_cutoff,
            match_score,
            mismatch_penalty,
            grid[h][0],
            grid[h][1],
        )
      
        if res:
            score = res[2]
            responses.append(res)
            #scores.append(score)
            hs.append(h)

            # exact match
            if score == 1.0:
                scores.append(score * len(res[1]))
                # cnt exact hit?
                #break
            else:
                scores.append(score)
        h += 1 

    if responses:
        idx = scores.index(max(scores))
        best_res = responses[idx]
        best_params = grid[hs[idx]]
        
        is_gapped_aln=False # to be removed

        candidate = best_res[0]
        
        gap_open_penalty, gap_extension_penalty = best_params[0], best_params[1]
        
        updated_reads = []        
        for read, aligner, ref_seq, ref_start in zip(best_res[1], best_res[5], best_res[3], best_res[4]):
            
            updated_read = update_read_info(
                                read,
                                candidate,
                                is_gapped_aln,
                                gap_open_penalty,
                                gap_extension_penalty,
                                aligner,
                                ref_seq,
                                ref_start,
                            )
                 
            updated_reads.append(updated_read)        
                     

        return candidate, updated_reads, gap_open_penalty, gap_extension_penalty
    else:
        return None                           
