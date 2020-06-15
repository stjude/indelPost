# cython: embedsignature=True

from .utilities import *
from pysam.libcfaidx cimport FastaFile
from pysam.libcbcf cimport VariantRecord, VariantFile


#cdef class Variant:
cdef class Variant:
    """This class accepts a VCF-style variant representation as input. 
    Equality holds between :class:`~indelpost.Variant` objects 
    if they are indentical after normalized. 

    Parameters
    ----------
    chrom : string
        contig name. 
    
    pos : integer
        1-based genomic position.
    
    ref : string
        VCF-style reference allele. 
    
    alt : string
        VCF-style alternative allele.
    
    reference : pysam.FastaFile
        reference FASTA file supplied as 
        `pysam.FastaFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.FastaFile>`__ object.
    
    Raises
    ------
    ValueError
        if the input locus is not defined in the reference or the input alleles contain letters 
        other than A/a, C/c, G/g, T/t, and N/n. 
        
    """
    cdef public str chrom, ref, alt
    cdef public int pos
    cdef readonly str _chrom
    cdef public FastaFile  reference
   
    #def __cinit__(self, str chrom, int pos, str ref, str alt, FastaFile reference):
    def __cinit__(self, chrom, pos, ref, alt, reference):
    #def __init__(self, chrom, pos, ref, alt, reference):
        self._chrom = str(chrom)
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.reference = reference

        self.chrom = self.__format_chrom_name(self._chrom, reference=reference)

        self.__validate()

    
    def __format_chrom_name(self, chrom, **kwargs):
        if kwargs.get("vcf", False):
            chrom_names = list(kwargs["vcf"].header.contigs)
        elif kwargs.get("reference", False):
            chrom_names = kwargs["reference"].references

        is_prefixed = True if chrom_names[0].startswith("chr") else False
        is_mt = True if "chrMT" in chrom_names or "MT" in chrom_names else False

        chrom = chrom.replace("chr", "")
        if chrom == "M" and is_mt:
            chrom = "MT"
        elif chrom == "MT" and not is_mt:
            chrom = "M"

        if is_prefixed:
            chrom = "chr" + chrom

        return chrom


    cpdef __validate(self):
        if not self.ref or not self.alt:
            raise ValueError("Allele may not be empty")

        if self.ref == self.alt:
            raise ValueError(
                "Not a variant: reference allele and alternate allele may not be identical"
            )

        bases = {"A", "C", "T", "G", "N", "a", "t", "c", "g", "n"}
        if not set(list(self.ref)) <= bases or not set(list(self.alt)) <= bases:
            raise ValueError("Allele contains char other than A/a, C/c, T/t, G/g, N/n")

        # check if contig is valid
        try:
            if not self.reference.fetch(self.chrom, self.pos - 1 , self.pos):
                raise ValueError("The locus is not defined in the reference")
        except:
            raise ValueError("The locus is not defined in the reference")


    @property
    def variant_type(self):
        """returns "I" if insertion, "D" if deletion, "S" if signle-nucleotide substitution, 
        and "M" if multi-nucleotide substitution. 
        """ 
        
        cdef int r_len, a_len
        cdef str var_type
            
        r_len, a_len = len(self.ref), len(self.alt)
        if r_len < a_len:
            var_type = "I"
        elif r_len > a_len:
            var_type = "D"
        elif r_len == a_len == 1:
            var_type = "S"
        else:
            var_type = "M"

        return var_type


    @property
    def is_del(self):
        """returns True if deletion
        """
        return self.variant_type == "D"


    @property
    def is_ins(self):
        """returns True if insertion
        """
        return self.variant_type == "I"

        
    @property
    def is_indel(self):
        """returns True if indel
        """
        return self.is_ins or self.is_del


    @property
    def indel_seq(self):
        """returns net inserted/deleted sequence. None if substitution.
        """
        if self.is_ins:
            return self.alt[len(self.ref) :]
        elif self.is_del:
            return self.ref[len(self.alt) :]
        else:
            return None


    def __eq__(self, other):
        cdef Variant i, j
        
        i, j = self.normalize(), other.normalize()

        equivalent = (
            i._chrom.replace("chr", "") == j._chrom.replace("chr", "")
            and i.pos == j.pos
            and j.ref.upper() == i.ref.upper()
            and i.alt.upper() == j.alt.upper()
        )

        return equivalent


    def __hash__(self):
        cdef Variant i = self.normalize() if self.is_indel else self
        hashable = (i._chrom, i.pos, i.ref, i.alt)

        return hash(hashable)


    def __dealloc__(self):
        pass

    @property
    def is_leftaligned(self):
        """returns True if left-aligned
        """
        if self.ref[-1].upper() != self.alt[-1].upper():
            return True

    
    @property
    def is_normalized(self):
        """returns True if left-aligned and parsimonious. 
        """
        if self.is_leftaligned:
            if len(self.ref) > 1 and len(self.alt) and (self.ref[0].upper() == self.alt[0].upper()):
                return False
            else:
                return True
        else:
            return False


    def normalize(self, inplace=False):
        """performes normalization on :class:`~indelpost.Variant` object.
        
        Parameters
        ----------

        inplace : bool
            returns None and normalizes this :class:`~indelpost.Variant` object if True. 
            Otherwise, returns a normalized copy of this object. Default to False.

        """
        cdef Variant i
        
        if inplace:
            i = self
        else:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)

        condition_1 = i.ref[-1].upper() == i.alt[-1].upper()
        while condition_1:
            left_base = i.reference.fetch(i._chrom, i.pos - 2, i.pos - 1)
            i.ref = left_base + i.ref[:-1]
            i.alt = left_base + i.alt[:-1]
            i.pos -= 1
            condition_1 = i.ref[-1].upper() == i.alt[-1].upper()

        condition_2 = i.ref[0].upper() == i.alt[0].upper()
        condition_3 = len(i.ref) > 1 and len(i.alt) > 1
        while condition_2 and condition_3:
            i.ref = i.ref[1:]
            i.alt = i.alt[1:]
            i.pos += 1
            condition_2 = i.ref[0].upper() == i.alt[0].upper()
            condition_3 = len(i.ref) > 1 and len(i.alt) > 1

        if inplace:
            return None
        else:
            return i
    
    
    def generate_equivalents(self):
        """generates non left-aligned copies of :class:`~indelpost.Variant` object.
        """ 
        
        i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference).normalize()
        pos, ref, alt = i.pos, i.ref, i.alt
        is_ins = i.is_ins
        
        res = [i]
        while self == i:
            right_base = i.right_flank(window=1)
            if is_ins:
                ref = alt[1]
                alt = alt[1:] + right_base
            else:
                alt = ref[1]
                ref = ref[1:] + right_base
            
            pos += 1 
            
            i = Variant(self.chrom, pos, ref, alt, self.reference)
            
            if self == i: 
                res.append(i)
        return res

    
    cpdef list query_vcf(self, VariantFile vcf, matchby="equivalence", window=50, as_dict=True):
        """finds VCF records matching this :class:`~indelpost.Variant` object.

        Parameters
        ----------
        vcf : pysam.VariantFile
            VCF file to be queried. 
            Supply as 
            `pysam.VariantFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantFile>`__ object.

        matchby : string
            "equivalence"
            "locus"
            "exact"
            
        window : integer
            VCF records  

        """

        cdef VariantRecord rec, hit

        matchbys = ["equivalence", "locus", "exact"]
        if not matchby in matchbys:
            raise ValueError("match by one of: %s" % matchbys)
        
        if self.variant_type == "S":
            leftaligned_pos, window = self.pos, 1
        else:
            leftaligned_pos = self.normalize().pos
        
        chrom = self.__format_chrom_name(self.chrom, vcf=vcf)
        
        try:
            records = to_flat_list(
                [
                    to_flat_vcf_records(rec)
                    for rec in vcf.fetch(chrom, leftaligned_pos - 1, leftaligned_pos - 1 + window)
                ]
            )
        except:
            return []

        hits = [
                record.orig
                for record in records
                if match_indels(Variant(self.chrom, record.pos, record.ref, record.alt, self.reference), self, matchby) ]
        
        if as_dict:
            hits = [
                {
                    "CHROM": hit.chrom,
                    "POS": hit.pos,
                    "ID": hit.id,
                    "REF": hit.ref,
                    "ALT": ",".join(list(hit.alts)),
                    "QUAL": hit.qual,
                    "FILTER": to_dict(hit.filter),
                    "INFO": to_dict(hit.info),
                    "FORMAT": to_dict(hit.format),
                    "SAMPLES": to_dict(hit.samples),
                }
                for hit in hits
            ]
        
        return hits


    def left_flank(self, window=50, normalize=False):
        if normalize:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)
        else:
            i = self

        lt_flank = i.reference.fetch(i.chrom, i.pos - window, i.pos)
        
        return lt_flank

    
    def right_flank(self, window=50, normalize=False):
        if normalize:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)
        else:
            i = self

        if i.variant_type == "I":
            rt_flank = i.reference.fetch(i.chrom, i.pos, i.pos + window)
        else:
            event_len = len(i.indel_seq) if "D" else len(i.ref)
            rt_flank = i.reference.fetch(i.chrom, i.pos + event_len, i.pos + event_len + window)

        return rt_flank

    
    def count_repeats(self, by_repeat_unit=True):
        if by_repeat_unit:
            seq = to_minimal_repeat_unit(self.indel_seq)
        else:
            seg = self.indel_seq

        lt_flank = self.left_flank()
        lr_repeat = repeat_counter(seq, lt_flank[::-1])
        rt_flank = self.right_flank()
        rt_repeat = repeat_counter(seq, rt_flank)
        
        return lr_repeat + rt_repeat
    

    def is_non_complex_indel(self):
        i = self.normalize()
        ref, alt = i.ref, i.alt
        if len(ref) == len(alt):
            return False
        
        if ref[0] != alt[0]:
            return False
        
        the_shorter = ref if i.is_ins else alt
        if len(the_shorter) > 1:
            return False
        
        return True
