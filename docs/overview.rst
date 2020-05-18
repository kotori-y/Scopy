..  -*- coding: utf-8 -*-

Overview
========
High-throughput screening (HTS) and virtual screening (VS) are widely applied in compounds screening and lead discovery. However, the frequent appearances of “noisy compounds” in the screening database, such as compounds with poor drug-likeness, compounds with poor selectivity or compounds with potential toxicity, have greatly weakened the efficiency of HTS and VS campaigns. The success of screening results critically depends on the quality of the available screening libraries. To construct a high-quality database, we developed Scopy (Screnning COmpounds in PYthon), an integrated negative design python library designed for screening out undesiable compounds in the early drug discovery. Scopy includes six modules, covering data preparation, screening filters, the calculation of scaffolds and descriptors, and the visualization analysis. The current version of Scopy can calculate 39 basic molecular properties, 3 comprehensive molecular evaluation scores, 2 types of molecular scaffolds, 6 types of substructure descriptors and 2 types of fingerprints. Screening rules such as drug-likeness rules (11 drug-likeness rules, 2 macro-cycle molecule rules and 2 building block rules), frequent hitter rules (4 assay interference substructure filters and 4 promiscuous compound substructure filters) and toxicophore filters (5 human related toxicity substructure filters, 3 environment related toxicity substructure filters and 3 comprehensive substructure filters) are provided in the Scopy library. Moreover, this library realized basic feature radar charts, feature-feature related scatter diagram, functional group marker gram and cloud gram four different visualization functions, which assists users in gaining a better understanding of the screening data. In conclusion, Scopy aims at providing an integrated analysis pipeline for molecule analysis, molecule screening, molecule optimization and model building. The Python package Scopy is designed by `CBDD Group`_ (Computational Biology & Drug Design Group), Xiangya School of Pharmaceutical Sciences, Central South University. 

.. _`CBDD Group`: http://home.scbdd.com/index.php?s=/Home/Index.html&t=english

Who uses Scopy?
~~~~~~~~~~~~~~~
Scopy is a comprehensive and uniform Python package which aims at providing an integrated pipeline for molecule analysis, molecule screening, molecule optimization and model building. Researchers who want to build a desirable screening database, evaluate initial hits quality, or need pertinent suggestions in molecular optimization, can use Scopy as a credible and efficient tool. We believe that, by the rational application of Scopy, users can draw useful knowledge and experience across many disciplines, thus decreasing the time and cost for drug research.
 
Motivation
~~~~~~~~~~
Scopy is intended to provide

-  Tools for pretreating molecules

-  Calculation of molecular physicochemical (PC) properties and screening drug-likeness filters
   
-  Detection of potential toxicophores
   
-  Detection of function groups related with frequent hitters

-  Calculation of substructure descriptors, fingerprints and scaffold

-  Visualization of molecular features and screening results
   
Feature overview
~~~~~~~~~~~~~~~~
Scopy contains six individual modules: (1) molecular preparation; (2) physicochemical properties and drug-likeness rules; (3) Molecular representations from substructure descriptors, fingerprints and molecular scaffolds; (4) FH substructure filters; (5) toxicophore filters and (6) the visualization analysis.

+-------------------------+----------------------------------------------------------------+-----------------------+
|Type                     |Detail                                                          |Description            |
+=========================+================================================================+=======================+
|Drug-likeness Filter     |Basic PC Properties                                             |                       |
|                         | - Molecular weight >>> MW                                      |                       |
|                         | - Molecular volume >>> Vol                                     |                       |
|                         | - Molecular density >>> Dense                                  |                       |
|                         | - Formal charge >>> fChar                                      |                       |
|                         | - Number of bonds >>> nBond                                    |                       |
|                         | - Number of atoms >>> nAtom                                    |                       |
|                         | - Number of heteroatoms >>> nHet                               |                       |
|                         | - Number of heavy atom >>> nHev                                |                       |
|                         | - Number of rotatable bonds >>> nRot                           |                       |
|                         | - Number of rigid bonds >>> nRig                               |                       |
|                         | - Molecular Flexibility >>> Flex                               |                       |
|                         | - Number of SSSR >>> nRing                                     |                       |
|                         | - logP >>> logP                                                |                       |
|                         | - logD >>> logD                                                |                       |
|                         | - logSw >>> logSw                                              |                       |
|                         | - Acid or Base >>> ab                                          |                       |
|                         | - pKa >>> pKa                                                  |                       |
|                         | - Molecular refraction >>>MR                                   |                       |
|                         | - Number of hydrogen bond donors >>> nHD                       |                       |
|                         | - Number of hydrogen bond acceptors >>> nHA                    |                       |
|                         | - Number of hydrogen bond donors& acceptors >>> nHB            |                       |
|                         | - Aromatic proportion >>> AP                                   |                       |
|                         | - Sp3 hybridized carbons/total carbon count >>> Fsp3           |                       |
|                         | - TPSA >>> TPSA                                                |                       |
|                         | - Size of biggest system ring >>> MaxRing                      |                       |
|                         | - Number of stero-centerss >>> nStero                          |                       |
|                         | - HetCarbonRatio >>> HetRatio                                  |                       |
|                         | - Number of single bonds >>> nSingle                           |                       |
|                         | - Number of double bonds >>> nDouble                           |                       |
|                         | - Number of triple bonds >>> nTriple                           |`ScoDruglikeness`_     |
|                         | - Number of Carbon atoms >>> nC                                |                       |
|                         | - Number of Boron atoms >>> nB                                 |                       |
|                         | - Number of Chlorine atoms >>> nCl                             |                       |
|                         | - Number of Bromine atoms >>> nBr                              |                       |
|                         | - Number of Iodine atoms >>> nI                                |                       |
|                         | - Number of Phosphor atoms >>> P                               |                       |
|                         | - Number of Sulfur atoms >>> nS                                |                       |
|                         | - Number of Oxygen atoms >>> nO                                |                       |
|                         | - Number of Nitrogen atoms >>> nN                              |                       |
|                         |Specific PC Properties                                          |                       |
|                         | - QED with average descriptor weights >>> QEDmean              |                       |
|                         | - QED with maximal descriptor weights >>> QEDmax               |                       |
|                         | - QED with using unit weights >>> QEDnone                      |                       |
|                         | - synthetic accessibility score >>> SAscore                    |                       |
|                         | - Natural product- likeness score >>> NPscore                  |                       |
+                         +----------------------------------------------------------------+                       +
|                         |Drug-likeness Rules (Small Molecule)                            |                       |
|                         | - Egan Rule                                                    |                       |
|                         | - Veber Rule                                                   |                       |
|                         | - LipinskiRule                                                 |                       |
|                         | - Pfizer Rule                                                  |                       |
|                         | - GSK Rule                                                     |                       |
|                         | - Oprea Rule                                                   |                       |
|                         | - Ghose Rule                                                   |                       |
|                         | - Xu Rule                                                      |                       |
|                         | - Ro4 Rule                                                     |                       |
|                         | - REOS Rule                                                    |                       |
|                         | - GoldenTriangle                                               |                       |
|                         |Drug-likeness Rules (Macro Molecule)                            |                       |
|                         | - BeyondRo5                                                    |                       |
|                         | - OralMacrocycles                                              |                       |
|                         |Drug-likeness Rules (Building Block)                            |                       | 
|                         | - Ro3 Rule                                                     |                       |
|                         | - Ro2 Rule                                                     |                       |
+-------------------------+----------------------------------------------------------------+-----------------------+
|Frequent Hitters Filter  |Assay interference                                              |                       |
|                         | - AlphaScreen_FHs(6)                                           |                       |
|                         | - Luciferase_Inhibitory(3)                                     |                       |
|                         | - Chelating(55)                                                |                       |
|                         | - Alarm_NMR(75)                                                |                       |
|                         | - Aggregator(311)                                              |                       |
|                         |Promiscuous compounds                                           |`ScoFH`_               |
|                         | - AlphaScreen_GST_FHs(34)                                      |                       |
|                         | - AlphaScreen_HIS_FHs(19)                                      |                       |
|                         | - PAINS(480)                                                   |                       |
|                         | - BMS(176)                                                     |                       |
+-------------------------+----------------------------------------------------------------+-----------------------+
|Toxicity Filter          |Human Toxicity                                                  |                       |
|                         | - Potential_Electrophilic(119)                                 |                       |
|                         | - LD50_oral(20)                                                |                       |
|                         | - Genotoxic_Carcinogenicity_Mutagenicity(117)                  |                       |
|                         | - NonGenotoxic_Carcinogenicity(23)                             |                       |
|                         | - Skin_Sensitization(155)                                      |                       |
|                         | - DNA_Binding(78)                                              |                       |
|                         |Comprehensive Toxicity                                          |`ScoTox`_              |
|                         | - NTD(105)                                                     |                       |
|                         | - SureChEMBL(165)                                              |                       |
|                         | - Toxicophores(154)                                            |                       |
|                         |Envrionment Toxicity                                            |                       |
|                         | - Acute_Aquatic_Toxicity(99)                                   |                       |
|                         | - Biodegradable(9)                                             |                       |
|                         | - NonBiodegradable(19)                                         |                       |
+-------------------------+----------------------------------------------------------------+-----------------------+
|Chemical Space Exploer   |Substurecture Descriptor                                        |                       |
|                         | - MACCS(167 bits)                                              |                       |
|                         | - **EFG(583 bits)**                                            |                       |
|                         | - **PubChem(881 bits)**                                        |                       |
|                         | - EState(79 bits)                                              |                       |
|                         | - **GhoseCrippen(110 bits)**                                   |                       |
|                         | - **IFG**                                                      |`ScoRepresent`_        |
|                         |Fingerptint                                                     |                       |
|                         | - Morgan(1024 *default*)                                       |                       |
|                         | - Daylight(2048 *default*)                                     |                       |
|                         |Framework                                                       |                       |
|                         | - Murcko Framework                                             |                       |
|                         | - Carbon Scaffold                                              |                       |
+-------------------------+----------------------------------------------------------------+-----------------------+

The Python programming language
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Python is a powerful programming language that allows simple and flexible representations of biochemical molecules, and clear and concise expressions of bioinformatics algorithms. Python has a vibrant and growing ecosystem of packages that Scopy uses to provide more features such as RDkit. In addition, Python is also an excellent “glue” language for putting together pieces of software from other languages which allows reuse of legacy code and engineering of high-performance algorithms. Equally important, Python is free, well-supported, and a joy to use. In order to make full use of Scopy, you will want to know how to write basic programs in Python. Among the many guides to Python, we recommend the documentation at https://www.python.org/

.. _`ScoDruglikeness`: ./modules/scopy.ScoDruglikeness.html
.. _`ScoFH`: ./modules/scopy.ScoFH.html
.. _`ScoTox`: ./modules/scopy.ScoTox.html
.. _`ScoRepresent`: ./modules/scopy.ScoRepresent.html