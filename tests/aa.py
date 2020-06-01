#!/usr/bin/env python3

import os
import pysam

from indeltools import Variant

tests = os.path.dirname(os.path.realpath(__file__))
reference = pysam.FastaFile(os.path.join(tests, "test.fa"))


from indeltools import Variant

indel_1 = Variant("chrZ", 1225, "C", "CAAGCCAAC", reference)
indel_2 = Variant("chrZ", 1228, "G", "GCCAACAAG", reference)

res = (indel_1 == indel_2)

print(res)
   
