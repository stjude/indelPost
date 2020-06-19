#from variant cimport Variant

#ctypedef fused str_or_list:
#    str
#    list

cdef tuple split(object data, str cigarstring, int target_pos, int string_pos, bint is_for_ref, bint reverse)

cdef tuple locate_indels(str cigarstring, int aln_start_pos)

#cdef tuple split_cigar(str cigarstring, int target_pos, int start)

#cdef tuple get_local_reference(Variant target, list pileup)

#cdef list cigar_lst  
#cdef int _size 

#cdef str cigar, event
#cdef int event_len, d_move, g_move
    
#cdef double [:] data_moves 
#cdef double [:] genome_moves 
    
#cdef int i, j    
