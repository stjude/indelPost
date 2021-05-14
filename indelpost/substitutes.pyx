#cython: embedsignature=True

from indelpost.variant cimport Variant, NullVariant
from pysam.libcalignmentfile cimport AlignmentFile

cdef class Substitute:
    """experimental class for substitute
    """
    def __cinit__(
        self, 
        Variant target,
        AlignmentFile bam,
        int window,
        bint exclude_duplicates.
        int downsample_threshold,
        int base_quality_threshold,
    ):
        
        self.target = target
        self.bam = bam
        self.window = window
        self.exclude_duplicates = exclude_duplicates


