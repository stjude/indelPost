# indelPost

[![Documentation Status](https://readthedocs.org/projects/indelpost/badge/?version=latest)](https://indelpost.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/indelpost.png)](https://badge.fury.io/py/indelpost)

indelPost is a Python library for indel processing via realignment and read-based phasing to resolve alignment ambiguities. By importing the library, 
users write their own scripts to solve alignment-sensitive problems such as:
* compare/integrate indels called differently by variant callers (e.g., complex indels)
* compare indel alignments in multiple mappings (e.g., match DNA indels to RNA-Seq for expression check)  
* construct a complex indel from a simple indel by read-based phasing    
* extract reads supporting the target indel from BAM file
* pull variant records matching the target indel from VCF file
* genotyping indels (the alleles may be aligned differently in the target BAM file)

Visit [documentation](https://indelpost.readthedocs.io/en/latest) for detail.

To install (require Linux with Python>=3.6 pre-installed):
```
pip install indelpost --no-binary indelpost --no-build-isolation
```

## Troubleshoot
If you get something like:
```
... may indicate binary incompatibility. Expected 88 from C header, got 72 from PyObject
```
or
```
AttributeError: module 'pysam.libcalignmentfile' has no attribute 'IteratorColumnAll'
```
try:
```
pip uninstall cython pysam indelpost
pip install cython pysam
pip install indelpost --no-binary indelpost --no-build-isolation
```

## Reference
Hagiwara K et al. (2022) indelPost: harmonizing ambiguities in simple and complex indel alignments. [Bioinformatics](https://doi.org/10.1093/bioinformatics/btab601)

## Contact
* kohei.hagiwara[AT]stjude.org 
