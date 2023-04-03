Limited support for SNV/MNVs
============================

For single- and multi-nucleotide variants (SNV/MNVs), :class:`~indelpost.Variant` object is fully supported.
:class:`~indelpost.VariantAlignment` object supports naive allele counting. 

.. image:: mnv.svg
   :width: 300
   :height: 30
   :align: center
|

"Naive" counting only includes aligned events ::

    import pysam
    from indelpost import Variant, VariantAlignment

    reference = pysam.FastaFile("/path/to/reference.fa")
    bam = pysam.AlignmentFile("/path/to/example.bam")

    v = Variant("chrN", 5, "GTC", "TAG", reference)
    valn = VariantAlignment(v, bam)

    cnt = valn.count_alleles()
    print(cnt)
    #(3, 2) the soft-clipped read (bottom) not included as target

