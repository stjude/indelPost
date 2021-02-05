# cython: embedsignature=True
# cython: profile=True

from .utilities import *
from .localn import make_aligner, align, findall_indels
from pysam.libcfaidx cimport FastaFile
from pysam.libcbcf cimport VariantFile


cdef class NullVariant:
    """This class represents a variant that does not exist. 
    Boolean expression evaluates to false. Alleles are always "NNN".
    
    Parameters
    ----------
    chrom : string
        contig name.

    pos : integer
        1-based genomic position

    reference : pysam.FastaFile
        reference FASTA file supplied as
        `pysam.FastaFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.FastaFile>`__ object.
    """
    def __cinit__(self, str chrom, int pos, FastaFile reference):
        self.chrom = chrom
        self.pos = pos
        self.ref = "NNN"
        self.alt = "NNN"
        self.reference = reference

    def __bool__(self):
        return False


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
    def __cinit__(self, str chrom, int pos, str ref, str alt, FastaFile reference):
        self._chrom = chrom
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
        ref_lst, alt_lst = list(self.ref), list(self.alt)
        if not set(ref_lst) <= bases or not set(alt_lst) <= bases:
            self.ref = "".join([base if base in bases else "N" for base in ref_lst])
            self.alt = "".join([base if base in bases else "N" for base in alt_lst])

        # check if contig is valid
        try:
            if not self.reference.fetch(self.chrom, self.pos - 1 , self.pos):
                raise ValueError("The locus is not defined in the reference")
        except:
            raise ValueError("The locus is not defined in the reference")


    @property
    def variant_type(self):
        """returns "I" if the net allele-length change is gain, 
        "D" if the net change is loss, 
        "S" if signle-nucleotide substitution (zero net change), 
        and "M" if multi-nucleotide substitution (zero net change). 
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
        """evaluates if :attr:`~indelpost.Variant.variant_type` is "D".
        """
        return self.variant_type == "D"


    @property
    def is_ins(self):
        """evaluates if :attr:`~indelpost.Variant.variant_type` is "I".
        """
        return self.variant_type == "I"

        
    @property
    def is_indel(self):
        """returns True if :attr:`~indelpost.Variant.variant_type` is "I" or "D".
        """
        return self.is_ins or self.is_del


    @property
    def indel_seq(self):
        """returns the inserted/deleted sequence for non-complex indels. None for substitutions.
        """
        if self.is_ins:
            return self.alt[len(self.ref) :]
        elif self.is_del:
            return self.ref[len(self.alt) :]
        else:
            return ""


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
        """returns True if the indel position is the smallest possible value (left-aligned). 
        """
        if self.ref[-1].upper() != self.alt[-1].upper():
            return True
        elif "N" in self.ref.upper() or "N" in self.alt.upper():
            return True

    
    @property
    def is_normalized(self):
        """returns True if left-aligned and the allele representations are minimal. 
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


    def _generate_equivalents_private(self):
        
        if self.is_non_complex_indel():
            return self.generate_equivalents()
        else:
            # define complex indel at the start and end of the deleted sequence 
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)
            j = Variant(self.chrom, self.pos + len(self.ref), self.ref, self.alt, self.reference)
            return [i, j]
    

    def query_vcf(self, VariantFile vcf, matchby="normalization", window=50, indel_only=True, as_dict=True):
        """finds VCF records matching this :class:`~indelpost.Variant` object.

        Parameters
        ----------
        vcf : pysam.VariantFile
            VCF file to be queried. 
            Supply as 
            `pysam.VariantFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantFile>`__ object.

        matchby : string
            "normalization" (default) matches by normalizing VCF records. "locus" matches by the normalized genomic locus.  
            "exact" finds the exact matche without normalization. 
            
        window : integer
            VCF records  

        """
        matchbys = ["normalization", "locus", "exact"]
        if not matchby in matchbys:
            raise ValueError("match by one of: %s" % matchbys)
        
        if self.variant_type == "S":
            leftaligned_pos, window = self.pos, 1
        else:
            leftaligned_pos = self.normalize().pos
        
        chrom = self.__format_chrom_name(self.chrom, vcf=vcf)
        
        searchable = vcf.fetch(chrom, leftaligned_pos - 1, leftaligned_pos - 1 + window)
        
        if not searchable:
            return []

        records = to_flat_list(
            [
                to_flat_vcf_records(rec)
                for rec in searchable
            ]
        )
        
        hits = [
                record.orig
                for record in records
                if match_indels(
                    Variant(self.chrom, record.pos, record.ref, record.alt, self.reference), 
                    self, 
                    matchby,
                    indel_only,
                )
        ]
        
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
        """extract the left-flanking reference sequence. See also :meth:`~indelpost.Variant.right_flank`.

        Parameters
        ----------
        window : int
            extract the reference sequence [variant_pos - window, variant_pos]. 
        normalize : bool
            if True, the normalized indel position is used as the end of the flanking sequence.
        """ 
        if normalize:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)
        else:
            i = self

        if i.is_non_complex_indel():
            pos = i.pos
        else:
            pos = i.pos - 1
        
        lt_flank = i.reference.fetch(i.chrom, max(0, pos - window), pos)
        
        return lt_flank

    
    def right_flank(self, window=50, normalize=False):
        """extract the right-flanking reference sequence. See also :meth:`~indelpost.Variant.left_flank`. 

        Parameters
        ----------
        window : int
            extract the reference sequence [variant_end_pos, variant_end_pos + window].
        normalize : bool
            if True, the normalized indel position is used as the start of the flanking sequence.
        """
        if normalize:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)
        else:
            i = self

        ref_lim = i.reference.get_reference_length(i.chrom)
        if i.is_non_complex_indel() and i.variant_type == "I":
            rt_flank = i.reference.fetch(i.chrom, i.pos, min(i.pos + window, ref_lim))
        else:
            if i.is_non_complex_indel() and i.variant_type == "D":
                event_len = len(i.indel_seq)
            else:
                event_len = len(i.ref) - 1  
            rt_flank = i.reference.fetch(i.chrom, i.pos + event_len, min(i.pos + event_len + window, ref_lim))

        return rt_flank

    
    def count_repeats(self, by_repeat_unit=True):
        """counts repeats in the flanking reference sequences. The search window is
        defined by :meth:`~indelpost.Variant.left_flank` and :meth:`~indelpost.Variant.right_flank`.
        
        Parameters
        ----------
        by_repeat_unit : bool
            count by the smallest tandem repeat unit. For example, the indel sequence "ATATATAT" has
            tandem units "ATAT" and "AT". The occurrent of the "AT" units is counted if True (default).
        """ 

        if self.is_non_complex_indel():
            seq = self.indel_seq
        else:
            seq = self.alt
        
        if by_repeat_unit:
            seq = to_minimal_repeat_unit(seq)

        lt_flank = self.left_flank()
        lr_repeat = repeat_counter(seq, lt_flank[::-1])
        rt_flank = self.right_flank()
        rt_repeat = repeat_counter(seq, rt_flank)
        
        return lr_repeat + rt_repeat
    

    def is_non_complex_indel(self):
        """returns True only if non-complex indel (False if complex indel or substitution).
        """
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


    def decompose_complex_variant(self, match_score=2, mismatch_penalty=2, gap_open_penalty=4, gap_extension_penalty=0):
        """returns a `list <https://docs.python.org/3/library/stdtypes.html#list>`__ of non-complex :class:`~indelpost.Variant` objects decomposed by the Smith-Waterman local alignment 
        with a given set of score/penalty.

        Parameters
        ----------
        match_score : int
            default to 2.
        mismatch_penalty : int
            default to 2.
        gap_open_penalty : int
            default to 4.
        gap_extension_penalty : int
            default to 0.    
        """
        if self.is_non_complex_indel():
            return [self]
        
        var = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference).normalize()
        
        lt_pos = var.pos - 1
        rt_pos = var.pos - 1 + len(var.ref)
        
        window = 100
        mut_seq = self.reference.fetch(var.chrom, lt_pos - window, lt_pos) + var.alt + self.reference.fetch(var.chrom, rt_pos, rt_pos + window)
        ref_seq = self.reference.fetch(var.chrom, lt_pos - window, lt_pos + window)

        aln = align(make_aligner(ref_seq, match_score, mismatch_penalty), mut_seq, gap_open_penalty, gap_extension_penalty) 
        
        genome_aln_pos = lt_pos + 1 - window + aln.reference_start
        
        indels, snvs = findall_indels(aln, genome_aln_pos, ref_seq, mut_seq, report_snvs=True)
        
        variants = []
        if indels:
            for idl in indels:
                padding_base = idl["lt_ref"][-1]
                if idl["indel_type"] == "D":
                    ref = padding_base + idl["del_seq"]
                    alt = padding_base
                else:
                    ref = padding_base
                    alt = padding_base + idl["indel_seq"]

                variants.append(Variant(self.chrom, idl["pos"], ref, alt, self.reference))
        
        if snvs:
            for snv in snvs:
                variants.append(Variant(self.chrom, snv["pos"], snv["ref"], snv["alt"], self.reference))
             
        return variants
