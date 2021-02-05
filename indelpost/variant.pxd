# cython: embedsignature=True

from pysam.libcfaidx cimport FastaFile
from pysam.libcbcf cimport VariantRecord, VariantFile

cdef class NullVariant:
    cdef public str chrom
    cdef readonly str ref, alt
    cdef public int pos
    cdef public FastaFile reference


cdef class Variant:
    cdef public str chrom, ref, alt
    cdef public int pos
    cdef public FastaFile reference
    cdef str _chrom
    cdef VariantFile vcf
    cdef VariantRecord rec, hit

    cpdef __validate(self)

