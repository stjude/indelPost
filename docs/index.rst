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
    
    pip install indelpost

or::
    
    git clone https://github.com/stjude/indelPost.git
    cd indelPost 
    python setup.py install

Installing indelPost also installs the dependency:

    - python>=3.6.0
    - numpy>=1.16.0
    - pysam>=0.15.0
    - cython>=0.29.12
    - ssw-py==0.2.6


.. note::
    
    indelPost is supported for Linux.

Contents
--------

.. toctree::
   :maxdepth: 2
   
   intro.rst
   api.rst
   examples.rst
   benchmark.rst


Indices and tables
------------------

Contetns:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
