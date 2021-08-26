Introduction
============

Small insertions/deletions (indels < 50 bp) can be represented differently between mappers and variant callers, or in the flanking sequence context. 
indelPost harmonizes these inter-mapper, inter-caller, and inter-read alignment ambiguities via realignment and read-based phasing to facilitate indel processing.

As an example, consider an inter-mapper ambiguity problem in DNA/RNA mappings.
Suppose that you want to test the expression of the DNA indels:

.. image:: dna.svg
   :width: 400
   :height: 400
   :align: center

|

Matching the DNA indels to the corresponding transcriptome data is non-trivial as the indels are aligned differently in the RNA-Seq:

.. image:: rna.svg
   :width: 400
   :height: 400
   :align: center

|

To test the expression, first prepare inputs::

    import pysam
    from indelpost import Variant, VariantAlignment
    
    # reference as pysam's FastaFile object
    reference = pysam.FastaFile("/path/to/reference.fa")
    
    # RNA-Seq BAM as pysam's AlignmentFile object
    bam = pysam.AlignmentFile("/path/to/RNASeq.bam")


Make a :class:`~indelpost.Variant` object to represent the target indel (the 2nd indel as example)::

    v = Variant("chrN", 9, "T", "TCCGTA", reference)

Next, make a :class:`~indelpost.VariantAlignment` object.
Upon creation, the :class:`~indelpost.VariantAlignment` object realigns the RNA-Seq reads to find the target indel::

    valn = VariantAlignment(v, bam)

Now, :meth:`~indelpost.VariantAlignment.count_alleles()` method counts the RNA-Seq reads supporing the target to check the expression::
    
    cnt = valn.count_alleles()
    print(cnt)  

The alignment equivalency of the first DNA indel, the second DNA indel, and the RNA indel can be tested::
 
    # first DNA indel
    first = Variant("chrN", 8, "G", "GCCAC", reference)
    valn_by_first = VariantAlignment(first, bam)

    # indel in RNA-Seq
    rna = Variant("chrN", 8, "G", "GCCACTCCGT", reference)
    valn_by_rna = VariantAlignment(rna, bam)
    
    # True
    valn == valn_by_first == valn_by_rna
   

These indels are equivalent because they have the same phased indel representation. :meth:`~indelpost.VariantAlignment.phase()` returns 
a :class:`~indelpost.Variant` representing the complex indel representation::
    
    v_phased = valn.phase() # same result for the other 2 objects.

    print(v_phased.chrom, v_phased.pos, v_phased.ref, v_phased.alt)
    
The last line outputs::
     
    "chrN", 9, "T", "CCACTCCGTA"

See more Usage  :ref:`examples` .

