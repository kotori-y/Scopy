.. Scopy documentation master file, created by
   sphinx-quickstart on Tue Nov 26 11:16:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Scopy's documentation
=========================
To decrease the four major problems existed in high-throughput screening (HTS), we have designed and develpoed the package Scopy(Screnning COmpounds in PYthon). Firstly, the poor drug-likeness of most hitters will increase workload for future work. Besides, the potential toxic compounds and frequent hitters (FH) would stain the screening result. Lastly, The chemical space of hitters may be narrow which may not involve the compound we intended.
Scopy supplied three filters: drug-likeness, toxicity and FH filter, but also a space analyser to explore chemical space of library and a visualizer to intuitively depict screening result. In drug-likeness filter, **43** physicochemical (PC) properties and **15** PC-driven rules, so called drug-likeness rules, were collected and implemented. In toxicity filter could screen **x** endpoints related to toxicity, and FH filter involves **x** endpoint. As to chemical space exploration, besides analyse framework, **7** fingerprins used to assess space from different angle.
The python package Scopy(Screnning COmpounds in PYthon) is designed by `CBDD Group`_ (Computational Biology & Drug Design Group), Xiangya School of Pharmaceutical Sciences, Central South University. 

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
