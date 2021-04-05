.. _Examples:

Examples
=========


Querying VCF file
-----------------
Below are *NF1* tumor suppressor gene indels regsitered in `COSMIC(v89) <https://cancer.sanger.ac.uk/cosmic>`__ . 
The VCF entries are juxtaposed with the alignments against GRCh38. 

.. image:: nf1.svg
   :width: 600
   :height: 50
   :align: center

|

To query the COSMIC VCF database for these indels::
    
    import pysam
    from indelpost import Variant
     
    reference = pysam.FastaFile("/path/to/GRCh38.fa")
    cosmic = pysam.VariantFile("/path/to/cosmic.v89.vcf(.gz)")

    # suppose your input is del2
    v = Variant("17", 31224665, "CC", "C", reference)
    
    # query by normalization (default)
    norm_hits = v.query_vcf(cosmic) 
    
    for hit in norm_hits:
        print(hit["INFO"]["CNT"]) #COSMIC counts for del1 and del2 

    # query for normalized genomic coordinate
    locus_hits = v.query_vcf(cosmic, matchby="locus")

    for hit in locus_hits:
        print(hit["INFO"]["CNT"]) #COSMIC counts for del1, del2, ins1, ins2, and in3
        
    # exact match (no normalization)
    exact_hit = v.query(cosmic, matchby="exact")
    
    print(exact_hit[0]["INFO"]["CNT"]) #COSMIC count for del2 only
    
    
    
Decomposing complex indels
--------------------------
Complex indels may be registered as a set of non-complex events in databases.
Decomposition may improve the searchability. Using `gnomAD(v3.0) <https://gnomad.broadinstitute.org>`__ as example:: 

    import pysam
    from indelpost import Variant

    reference = pysam.FastaFile("/path/to/GRCh38.fa")
    gnomad = pysam.VariantFile("/path/to/gnomad.3.0.vcf(.gz))

    v = Variant("chr1", 114299169, "CAGTGA", "TCTCT", reference)

    # no hits -> empty list
    no_hits = v.query_vcf(gnomad)

    # list of non-complex Variant objects 
    decomposed = v.decompose_complex_variant()
    
    for d in decomposed:
        
        # chr1 114299168 A AT
        # chr1 114299169 CAG C
        # chr1 114299173 G C
        # chr1 114299174 A T
        
        print(d.chrom, d.pos, d.ref, d.alt)
      
        if d.query_vcf(gnomad):
        
            # do something...
            
    








 
