from indelpost.variant cimport Variant

cdef class Contig:
    cdef Variant target
    cdef list pileup, targetpileup
    cdef public bint qc_passed, failed
    cdef int donwsample_lim
    cdef object lt_genomic_index, rt_genomic_index
    cdef int start, end
    cdef str lt_reference_seq, rt_reference_seq
    cdef double low_consensus_thresh
    cdef public int is_target_right_aligned
    cdef public object contig_dict
    cdef public str lt_consensus_seq, rt_consensus_seq, indel_seq
    cdef public list lt_consensus_scores, rt_consensus_scores, mismatches, non_target_indels, gaps
    cdef public int mapq
    cdef public double low_qual_mapping_rate
    cdef public dict qc_stats
