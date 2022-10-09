.. indelPost documentation master file, created by
   sphinx-quickstart on Tue Jun  2 18:31:25 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

indelPost: Python library for simple and complex indels. 
===================================================================

:Author: Kohei Hagiwara
:Date: |today|
:Version: |version|


To install::
    
    git clone https://github.com/stjude/indelPost.git
    cd indelPost 
    python setup.py install

Installing indelPost also installs the dependency:

    - numpy>=1.16.0
    - pysam>=0.15.0
    - cython>=0.29.12


.. note::
    
    indelPost is supported for Linux with Python>=3.8.

Contents
--------

.. toctree::
   :maxdepth: 2
   
   intro.rst
   api.rst
   examples.rst
   benchmark.rst

* :ref:`genindex`
