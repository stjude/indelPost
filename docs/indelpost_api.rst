Introduction
============

Small indels (< 51 bp) are frequent mutations and their representations may not be unique in short read NGS data, which complicates downstream analyses.
indelPost is a library to facilitate the processing of indel calls/alignments by resolving the alignment ambiguities 
via realignment and local-phasing.

As a quick example, suppose that you want to know if the following DNA indels are expressed:

.. image:: dna.svg
   :width: 400
   :height: 400
   :align: center

|

But, the indels are aligned differently in the RNA-Seq:

.. image:: rna.svg
   :width: 400
   :height: 400
   :align: center

|

To test the expression, write a script like this::

    import pysam
    from indelpost import Variant, VariantAlignment
    
    # reference as pysam's FastaFile object
    reference = pysam.FastaFile("/path/to/reference.fa")
    
    # RNA-Seq BAM as pysam's AlignmentFile object
    bam = pysam.AlignmentFile("/path/to/RNASeq.bam")

    # create indelpost's Variant object (first DNA indel used here)
    v1 = Variant("chrN", 8, "G", "GCCAC", reference)

    # RNA-Seq pileup is realigned to find the input indel
    aln1 = VariantAlignment(v1, bam)

    # check the allele count in the returned tuple (#REF allele, #ALT allele)
    aln1.count_alleles()

Check the identity of these seemingly different alignments::
 
    # second DNA indel
    v2 = Variant("chrN", 9, "T", "TCCGTA", reference)
    aln2 = VariantAlignment(v2, bam)

    
    # indel in RNA-Seq
    v3 = Variant("chrN", 8, "G", "GCCACTCCGT", reference)
    aln3 = VariantAlignment(v3, bam)

    
    # True
    are_identical = (aln1 == aln2 == aln3)
    
:class:`~indelpost.VariantAlignment` objects are indentical if :meth:`~indelpost.VariantAlignment.phase` 
returns the same :class:`~indelpost.Variant` object::

    phased1 = aln1.phase() 
    phased2 = aln2.phase()
    phased3 = aln3.phase()

    phased1 == phased2 == phased3 == Variant("chrN", 9, "T", "CCACTCCGTA", reference)

See more Usage  :ref:`examples` .

API
===

indel processing on the reference
---------------------------------
.. autoclass:: indelpost.Variant()
   :members:
   :exclude-members: chrom, pos, ref, alt, reference

.. autoclass:: indelpost.NullVariant()

indel processing on alignment files (BAM)
-----------------------------------------
.. autoclass:: indelpost.VariantAlignment()
   :members:

.. autoclass:: indelpost.Contig()
   :members:
   :exclude-members: contig_dict, failed, gaps, indel_seq, is_target_right_aligned, 
                     low_qual_mapping_rate, lt_consensus_scores, lt_consensus_seq,
                     lt_target_block_consensus_scores, lt_target_block_consensus_seq,
                     mapq, mismatches, non_target_indels, qc_passed, qc_stats, rt_consensus_scores,
                     rt_consensus_scores, rt_consensus_seq, rt_target_block_consensus_seq, rt_target_block_consensus_scores,
                     splice_pattern

.. autoclass:: indelpost.FailedContig()
   :members:
   :exclude-members: target_not_found, is_low_quality, failed_anyway
