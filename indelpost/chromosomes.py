#!/usr/bin/env python3
import pysam

CANONICALS = [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]


def chromsome_info(filename, preset):
    """Get chromosome name list
    
    Args:
        filename (file): file to be checked for chromosome name format
        format (str): specify input file format: fasta, bam, vcf
    Returns:
        chrom_names (lst): a list of chromosome names 
    Raises:
        Exception if invalid format is provided
    """
    if preset == "fasta":
        genome = pysam.FastaFile(filename)
        chrom_names = genome.references
    elif preset == "bam":
        bam_header_dict = pysam.AlignmentFile(filename).header["SQ"]
        chrom_names = [d["SN"] for d in bam_header_dict]
    elif preset == "vcf":
        vcf = pysam.VariantFile(filename)
        chrom_names = list(vcf.header.contigs)
    else:
        raise Exception("Valid formats are fasta, bam, or vcf")

    return chrom_names


def chrom_format(chrom, filename, preset="fasta"):
    is_prefixed, is_mt = chromosome_name_format(filename, preset=preset)
    return format_chromosome_name(chrom, is_prefixed=is_prefixed, is_mt=is_mt)


def chromosome_name_format(filename, preset="fasta"):
    """Check chromosome name format
    
    Args:
        filename (file): file to be checked for chromosome name format
        format (str): specify input file format: fasta, bam, vcf. Default to fasta
    Returns:
        tuple (is_prefixed, is_mt)
            is_prefixed (bool): True if chromosome names are chr-prefixed
            is_mt (bool): True if mitochondoria chromosome name is chrMT or MT
    """
    chrom_names = chromsome_info(filename, preset)

    is_prefixed = True if chrom_names[0].startswith("chr") else False
    is_mt = True if "chrMT" in chrom_names or "MT" in chrom_names else False

    return is_prefixed, is_mt


def format_chromosome_name(chrom, is_prefixed=True, is_mt=True):
    """Convert chr name to specified format

    Args:
        chrom (str): chromosome name
        is_prefixed (bool): True if to be chr-prefixed. Default to True
        is_mt (bool): True if to be chrMT or MT. Default to True
    Returns:
        chrom (str): formatted name
    """
    chrom = chrom.replace("chr", "")
    if chrom == "M" and is_mt:
        chrom = "MT"
    elif chrom == "MT" and not is_mt:
        chrom = "M"

    if is_prefixed:
        chrom = "chr" + chrom

    return chrom


def is_canonicalom(chr):
    """Ask if chromosome is 1-22, X, Y, or MT
    
    Args:
        chrom (str): chromosome name
    Returns:
        bool: True if 1-22, X, Y, or MT
    """
    chrom = chrom.replace("chr", "")

    return chrom in CANONICALS


def chromosome_order(filename, preset="fasta"):
    """Store chromosome order

    Args:
        filename (file): file to be checked for chromosome name format
        format (str): specify input file format: fasta, bam, vcf. Default to fasta
    Returns:
        dict: {chrosomome name (str)  : order (int)}
    """
    chrom_names = chromsome_info(filename, preset)
    return {chrom: i for i, chrom in enumerate(chrom_names)}


if __name__ == "__main__":
    fasta = "/rgs01/resgen/prod/tartan/index/reference/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa"
    a = chromosome_order(fasta, preset="fasta")
    print(a["chr3"])
