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

To install (require Linux with Python>=3.6 pre-installed):
```
pip install indelpost --no-binary indelpost --no-build-isolation
```

Visit [documentation](https://indelpost.readthedocs.io/en/latest) for detail.

