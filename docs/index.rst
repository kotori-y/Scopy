.. Scopy documentation master file, created by
   sphinx-quickstart on Tue Nov 26 11:16:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Scopy's documentation
===========================
High-throughput screening (HTS) and virtual screening (VS) are widely applied in compounds screening and lead discovery. However, the frequent appearances of “noisy compounds” in the screening database, such as compounds with poor drug-likeness, compounds with poor selectivity or compounds with potential toxicity, have greatly weakened the efficiency of HTS and VS campaigns. The success of screening results critically depends on the quality of the available screening libraries. To construct a high-quality database, we developed Scopy (Screnning COmpounds in PYthon), an integrated negative design python library designed for screening out undesiable compounds in the early drug discovery. Scopy includes six modules, covering data preparation, screening filters, the calculation of scaffolds and descriptors, and the visualization analysis. The current version of Scopy can calculate 39 basic molecular properties, 3 comprehensive molecular evaluation scores, 2 types of molecular scaffolds, 6 types of substructure descriptors and 2 types of fingerprints. Screening rules such as drug-likeness rules (11 drug-likeness rules, 2 macro-cycle molecule rules and 2 building block rules), frequent hitter rules (4 assay interference substructure filters and 4 promiscuous compound substructure filters) and toxicophore filters (5 human related toxicity substructure filters, 3 environment related toxicity substructure filters and 3 comprehensive substructure filters) are provided in the Scopy library. Moreover, this library realized basic feature radar charts, feature-feature related scatter diagram, functional group marker gram and cloud gram four different visualization functions, which assists users in gaining a better understanding of the screening data. In conclusion, Scopy aims at providing an integrated analysis pipeline for molecule analysis, molecule screening, molecule optimization and model building. The Python package Scopy is designed by `CBDD Group`_ (Computational Biology & Drug Design Group), Xiangya School of Pharmaceutical Sciences, Central South University. 

.. _`CBDD Group`: http://home.scbdd.com/index.php?s=/Home/Index.html&t=english

.. toctree::
   :maxdepth: 3

   overview
   user_guide
   application

.. toctree::
   :maxdepth: 3

   modules
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. figure:: /image/logocbdd.png
    :width: 500px
    :align: center