# cython: embedsignature=True

from pysam.libcfaidx cimport FastaFile
from pysam.libcbcf cimport VariantRecord, VariantFile


cdef class Variant:
    cdef public str chrom, ref, alt
    cdef public int pos
    cdef readonly str _chrom
    cdef public FastaFile  reference
    cdef VariantFile vcf
    cdef VariantRecord rec, hit

    cpdef __validate(self)
    
    #cpdef list query_vcf(self, VariantFile vcf, str matchby=*, int window=*, bint as_dict=*)
