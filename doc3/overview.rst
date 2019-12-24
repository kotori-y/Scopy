..  -*- coding: utf-8 -*-

Overview
========
To decrease the four major problems existed in high-throughput screening (HTS), we have designed and develpoed the package Scopy(Screnning COmpounds in PYthon). Firstly, the poor drug-likeness of most hitters will increase workload for future work. Besides, the potential toxic compounds and frequent hitters (FH) would stain the screening result. Lastly, The chemical space of hitters may be narrow which may not involve the compound we intended.
Scopy supplied three filters: drug-likeness, toxicity and FH filter, but also a space analyser to explore chemical space of library and a visualizer to intuitively depict screening result. In drug-likeness filter, **43** physicochemical (PC) properties and **15** PC-driven rules, so called drug-likeness rules, were collected and implemented. In toxicity filter could screen **x** endpoints related to toxicity, and FH filter involves **x** endpoint. As to chemical space exploration, besides analyse framework, **7** fingerprins used to assess space from different angle.
The python package Scopy(Screnning COmpounds in PYthon) is designed by `CBDD Group`_ (Computational Biology & Drug Design Group), Xiangya School of Pharmaceutical Sciences, Central South University. 

.. _`CBDD Group`: http://home.scbdd.com/index.php?s=/Home/Index.html&t=english

Who uses Scopy?
~~~~~~~~~~~~~~~

For those researchers from different biomedical fields, the Scopy package can be used to enhance the chance of finding compounds with an acceptable ADMET profile. Scopy will be helpful when finding the method to reduce the space of compound space.

 
Motivation
~~~~~~~~~~
Scopy is intended to provide

-  Analyse the physicochemical (PC) properties and filter compounds based on PC-derived rules
   
-  Search for the presence of toxicophores
   
-  Flag unwanted reactive chemical groups

-  Compute fingerprints retrieved from fragments

-  Provide an intuitive visualization
   
Feature overview
~~~~~~~~~~~~~~~~

Scopy cantain five parts: **drug-likeness filter** to analyse the physicochemical (PC) properties and filter compounds based on PC-derived rules, **structure alert filter** to search for the presence of toxicophores and flag unwanted reactive chemical groups, **fingerprints calulater** to compute fingerprints retrieved from fragments, **visualization tool** to provide an intuitive visualization, **pretreat** to solve some molecules can not be recognize. The table beblow shows each the detail of each part. Roughly, there are 41 PC properties and 15 PC-derived rules to be analysed. Besides, Scopy also contains 2,167 unexpected groups  covering 21 endpoints. In addtion, 7 fingerprints retrieved from fragments are supplied.

+-------------------------+-----------------------------------------------------------------------------------------------+-----------------------+
|Type                     |Detail                                                                                         |Description            |
+=========================+===============================================================================================+=======================+
|Drug-likeness filter     |PC Properties                                                                                  |                       |
|                         | - Molcular Weight >>> MW                                                                      |                       |
|                         | - Molcular Volume >>> MolVol                                                                  |                       |
|                         | - Molcular Density >>> Dense                                                                  |                       |
|                         | - Formal Charge >>> fChar                                                                     |                       |
|                         | - Number of bonds >>> nBond                                                                   |                       |
|                         | - Number of atoms >>> nAtom                                                                   |                       |
|                         | - Number of heteroatoms >>> nHet                                                              |                       |
|                         | - Number of heavy atom >>> nHev                                                               |                       |
|                         | - Number of rotable bonds >>> nRot                                                            |                       |
|                         | - Number of rigid bonds >>> nRig                                                              |                       |
|                         | - Number of SSSR >>> nRing                                                                    |                       |
|                         | - logP >>> logP                                                                               |                       |
|                         | - logD >>> logD                                                                               |                       |
|                         | - logSw >>> logSw                                                                             |                       |
|                         | - Acid or Base >>> ab                                                                         |                       |
|                         | - pKa >>> pKa                                                                                 |                       |
|                         | - QED with average descriptor weights >>> QEDmean                                             |                       |
|                         | - QED with maximal descriptor weights >>> QEDmax                                              |                       |
|                         | - QED with using unit weights >>> QEDnone                                                     |                       |
|                         | - Molecular refraction >>>MR                                                                  |                       |
|                         | - Number of hydrogen bond donors >>> nHD                                                      |                       |
|                         | - Number of hydrogen bond acceptors >>> nHA                                                   |                       |
|                         | - Number of hydrogen bond donors& acceptors >>> nHB                                           |                       |
|                         | - Aromatic proportion >>> AP                                                                  |                       |
|                         | - Sp3 hybridized carbons/total carbon count >>> Fsp3                                          |                       |
|                         | - TPSA >>> TPSA                                                                               |                       |
|                         | - Number of atoms involved in the biggest system ring >>> MaxRing                             |                       |
|                         | - Number of Sterocenterss >>> nStero                                                          |`druglikeness`_        |
|                         | - HetCarbonRatio >>> HetRatio                                                                 |                       |
|                         | - synthetic accessibility score >>> SAscore                                                   |                       |
|                         | - Natural product- likeness score >>> NPscore                                                 |                       |
|                         | - Number of single bonds >>> nSingle                                                          |                       |
|                         | - Number of double bobds >>> nDouble                                                          |                       |
|                         | - Number of triple bonds >>> nTriple                                                          |                       |
|                         | - Number of Carbon atoms >>> nC                                                               |                       |
|                         | - Number of Boron atoms >>> nB                                                                |                       |
|                         | - Number of Chlorin atoms >>> nCl                                                             |                       |
|                         | - Number of Bromine atoms >>> nBr                                                             |                       |
|                         | - Number of Iodine atoms >>> nI                                                               |                       |
|                         | - Number of Phosphor atoms >>> P                                                              |                       |
|                         | - Number of Sulfur atoms >>> nS                                                               |                       |
|                         | - Number of Oxygen atoms >>> nO                                                               |                       |
|                         | - Number of Nitrogen atoms >>> nN                                                             |                       |
+                         +-----------------------------------------------------------------------------------------------+                       +
|                         |PC- derived RRules                                                                             |                       |
|                         | - Egan Rule     0<=tPSA<=132; -1<=logP<=6                                                     |                       |
|                         | - Veber Rule    nRot<= 10; TPSA<=140; nHB<=12                                                 |                       |
|                         | - LipinskiRule  MW<=500; logP<=5, nHD<=5, nHA<=10                                             |                       |
|                         | - BeyondRo5   MW<=1000; -2<=logP<=10; nHD<=6, nHA<=15; tPSA<=250; nRot<=20                    |                       |
|                         | - Pfizer Rule     logP>3; TPSA<75                                                             |                       |
|                         | - GSK Rule        MW<=400; logP<=4                                                            |                       |
|                         | - OralMacrocycles     MW<1000; logP<10; nHD<5; PSA<250                                        |                       |
|                         | - Oprea Rule   nRing>=3,nRig>=18,nRot>=6                                                      |                       |
|                         | - Ghose Rule      -0.4<logP<5.6; 160<MW<480; 40<MR<130; 20<nAtom<70                           |                       |
|                         | - Xu Rule     nHD<=5; nHA<=10; 3<=nRot<= 35; 1<=nring<=7; 10<=nhev<=50                        |                       |
|                         | - Ro4 Rule    MW<=400; logP<=4; nHDv=4; NHA<=8; PSAv=120                                      |                       |
|                         | - Ro3 Rule    MW<=300; -3<=logP<=3; nHD<=3; nHA<=6; PSA<=60                                   |                       |
|                         | - Ro2 Rule    MWv=200; logP<=2; nHD<=2; nHA<=4                                                |                       |
|                         | - REOS Rule    200<=MW<=500; -5<=logP<=5; nHD<=5; nHA<=10; nRot<=8; TPSA<=150; -4<=fChar<=4   |                       |
|                         | - GoldenTriangle 200<=MW<=500; -2<=logD<=5                                                    |                       |
+-------------------------+-----------------------------------------------------------------------------------------------+-----------------------+
|Structure alert filter   |Unwanted fragments                                                                             |                       |
|                         | - Acute_Aquatic_Toxicity(99)                                                                  |                       |
|                         | - AlphaScreen_FHs(6)                                                                          |                       |
|                         | - AlphaScreen_GST_FHs(34)                                                                     |                       |
|                         | - AlphaScreen_HIS_FHs(19)                                                                     |                       |
|                         | - Biodegradable(9)                                                                            |                       |
|                         | - BMS(180)                                                                                    |                       |
|                         | - Chelating(55)                                                                               |                       |
|                         | - Developmental_Mitochondrial(12)                                                             |                       |
|                         | - Genotoxic_Carcinogenicity_Mutagenicity(117)                                                 |                       |
|                         | - Idiosyncratic(35)                                                                           |                       |
|                         | - LD50_oral(20)                                                                               |                       |
|                         | - Luciferase_Inhibitory(3)                                                                    |                       |
|                         | - NonBiodegradable(19)                                                                        |`structure_alert`_     |
|                         | - NonGenotoxic_Carcinogenicity(23)                                                            |                       |
|                         | - NTD(105)                                                                                    |                       |
|                         | - Pains(480)                                                                                  |                       |
|                         | - Potential_Electrophilic(119)                                                                |                       |
|                         | - Promiscuity(177)                                                                            |                       |
|                         | - Reactive_Unstable_Toxic(335)                                                                |                       |
|                         | - Skin_Sensitization(155)                                                                     |                       |
|                         | - SureChEMBL(165)                                                                             |                       |
|                         | - BMS(176)                                                                                    |                       |
|                         | - NTD(105)                                                                                    |                       |
|                         | - Alarm_NMR(75                                                                                |                       |
|                         | - Frequent_Hitters(15)                                                                        |                       |
|                         | - Aggregators(311)                                                                            |                       |
|                         | - Toxicophores(154)                                                                           |                       |
+-------------------------+-----------------------------------------------------------------------------------------------+-----------------------+
|Fingerprints Calculater  | - MACCS(167- bits)                                                                            |                       |
|                         | - Morgan(2^n- bits,- 1024- set- as- default)                                                  |                       |
|                         | - **EFG(583- bits)**                                                                          |                       |
|                         | - Daylight(2^n- bits,- 2048- set- as- default)                                                |`fingerprints`_        |
|                         | - **PubChem(881- bits)**                                                                      |                       |
|                         | - EState(79- bits)                                                                            |                       |
|                         | - **GhoseCrippen(110- bits)**                                                                 |                       |
+-------------------------+-----------------------------------------------------------------------------------------------+-----------------------+




The Python programming language
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python is a powerful programming language that allows simple and flexible representations of biochemical molecules, and clear and concise expressions of bioinformatics algorithms. Python has a vibrant and growing ecosystem of packages that Scopy uses to provide more features such as RDkit. In addition, Python is also an excellent “glue” language for putting together pieces of software from other languages which allows reuse of legacy code and engineering of high-performance algorithms. Equally important, Python is free, well-supported, and a joy to use. In order to make full use of Scopy, you will want to know how to write basic programs in Python. Among the many guides to Python, we recommend the documentation at https://www.python.org/