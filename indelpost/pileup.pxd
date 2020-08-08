from indelpost.variant cimport Variant
from pysam.libcalignmentfile cimport AlignmentFile

cdef tuple make_pileup(
    Variant target, 
    AlignmentFile bam, 
    bint exclude_duplicates, 
    int window, 
    int downsamplethresh, 
    int basequalthresh
)
