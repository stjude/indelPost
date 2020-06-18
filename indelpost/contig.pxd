from variant cimport Variant

cdef class Contig:
    cdef Variant target
    cdef list pileup, targetpileup
    cdef bint qc_passed, failed
    cdef int donwsample_lim, mapq
    cdef dict lt_genomic_index, rt_genomic_index, genome_indexed_contig
    cdef int start, end
    cdef str lt_reference_seq, rt_reference_seq
    cdef str lt_consensus_seq, rt_consensus_seq
    cdef list lt_consensus_scores, rt_consensus_scores
    cdef str indel_seq
    cdef list non_target_indels, mismatches, gaps
