
#ctypedef fused str_or_list:
#    str
#    list

cdef tuple split(object data, str cigarstring, int target_pos, int string_pos, bint is_for_ref, bint reverse)
#cpdef tuple qsplit(list data, str cigarstring, int target_pos, int string_pos, bint is_for_ref, bint reverse)

#cdef list cigar_lst  
#cdef int _size 

#cdef str cigar, event
#cdef int event_len, d_move, g_move
    
#cdef double [:] data_moves 
#cdef double [:] genome_moves 
    
#cdef int i, j    
