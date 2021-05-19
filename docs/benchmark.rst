Benchmarking
===============

:class:`~indelpost.VariantAlignment` object construction task for indels called by `GATK4 HaplotypeCaller(4.0.2.1) <https://gatk.broadinstitute.org/hc/en-us/articles/360037059732-HaplotypeCaller>`__ 
on 75 whole-exome sequencing data `SJC-DS-1003 <https://platform.stjude.cloud/data/cohorts>`__. `pandas DataFrame <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html>`__ used to store the indel calls. 

.. image:: resource_usage.png
   :width: 600
   :align: center


.. tip::
    
    Consider split data into smaller chunks and run in parallel for genome-wide indel datasets.

