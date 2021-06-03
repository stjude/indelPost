# cython: embedsignature=True

#from pysam.libcbcf cimport VariantRecord, VariantFile
from pysam.libcalignmentfile cimport AlignmentFile
from indelpost.variant cimport Variant
from indelpost.contig cimport Contig

cdef class VariantAlignment:
    cdef Variant target, __target
    cdef readonly AlignmentFile bam
    cdef int window, retarget_window, mapqthresh
    cdef int downsamplethresh, basequalthresh, match_score, _observed_pos
    cdef int mismatch_penalty, gap_open_penalty, gap_extension_penalty
    cdef float retarget_cutoff, __sample_factor,
    cdef bint exclude_duplicates, auto_adjust_extension_penalty
    cdef list __pileup
    cdef readonly is_spurious_overhang
    cdef readonly Contig contig

    cdef __parse_pileup(self, Contig contig=*, bint retargeted=*)
