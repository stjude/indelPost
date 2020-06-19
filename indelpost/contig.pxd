from variant cimport Variant

cdef class Contig:
    cdef Variant target
    cdef list pileup, targetpileup
    cdef public bint qc_passed, failed
    cdef int donwsample_lim, mapq
    cdef object lt_genomic_index, rt_genomic_index, genome_indexed_contig
    cdef int start, end
    cdef str lt_reference_seq, rt_reference_seq
    cdef public str lt_consensus_seq, rt_consensus_seq
    cdef public list lt_consensus_scores, rt_consensus_scores
    cdef public str indel_seq
    cdef list non_target_indels, mismatches, gaps
