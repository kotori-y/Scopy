.. Scopy documentation master file, created by
   sphinx-quickstart on Tue Nov 26 11:16:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Scopy's documentation
=========================
The python package Scopy(Screnning COmpounds in PYthon) is designed by `CBDD Group`_ (Computational Biology & Drug Design Group), Xiangya School of Pharmaceutical Sciences, Central South University. The noise, existed in compound library, may interfere the screening result, so called False Positive. To help select quality hit compounds, a filter is necessary before screenig and(or) purchasing. Scopy is a such filter that can not only analyse the physicochemical (PC) properties and filter compounds based on PC-derived rules, but search for the presence of toxicophores (potentially toxic chemical groups). Scopy can also flag unwanted reactive chemical groups, e.g. Pan Assay Interference Compounds (PAINS) or aggregators Our aims at enhancing the chances of finding compounds with an acceptable absorption, distribution, metabolism, excretion and toxicity (ADMET) profile. Roughly, there are 41 PC properties and 15 PC-derived rules to be analysed. Besides, Scopy also contains 2,167 unexpected groups  covering 21 endpoints. In addtion, 7 fingerprints retrieved from fragments are supplied.

.. _`CBDD Group`: http://home.scbdd.com/index.php?s=/Home/Index.html&t=english

.. toctree::
   :maxdepth: 3

   overview
   user_guide


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
