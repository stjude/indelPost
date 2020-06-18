#!/usr/bin/env python3
import time
import pysam
import indelpost as ip


fasta = "/research/rgs01/resgen/prod/tartan/runs/ad_hoc/reference_import-Tm4OgjZi/output/reference/Homo_sapiens/GRCh37-lite/FASTA/GRCh37-lite.fa"
#fasta = "/rgs01/resgen/prod/tartan/index/reference/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa"
reference = pysam.FastaFile(fasta)


bam = "/rgs01/resgen/prod/tartan/index/data/SJLIFE/DeepSeqhg19/SJOS042848_G1/VALCAP/bam/SJOS042848_G1.bam"
bam = pysam.AlignmentFile(bam)


vcf = "/home/khagiwar/git/OutLier_orig/data/dbsnp_v151_grch37.vcf.gz"
vcf = pysam.VariantFile(vcf)


chrom = "2"
chrom = "20"
pos = 25457161
pos = 57415542

ref = "AAAT"
ref = "CGAGCCTGAGACCGCCCCCACCACT"
alt = "A"
alt = "C"


start = time.time()
var = ip.Variant(chrom, pos, ref, alt, reference)

a = var.query_vcf(vcf, matchby="locus")

aln = ip.VariantAlignment(var, bam, exclude_duplicates=False)
b = aln.count_alleles()
print(b)
print(time.time() - start)
