..  -*- coding: utf-8 -*-

Getting Started with Scopy
==========================
This document intends to provide users with the basic operation methods of Scopy. If you find any mistake or have suggestions for improvements, please either fix them in the source document (the .py file) or send to the mailing list: oriental-cds@163.com and kotori@cbdd.me.

Installing the Scopy package
-----------------------------
Scopy has been successfully tested on Linux, OSX and Windows systems under Python3.6 and Python3.7

Dependencies
~~~~~~~~~~~~
.. code-block:: python

	rdkit
	numpy
	matplotlib

Install with source
~~~~~~~~~~~~~~~~~~~
>>> git clone git@github.com:kotori-y/Scopy.git && cd scopy
>>> [sudo] python setup.py install

Install with conda
~~~~~~~~~~~~~~~~~~~
>>> conda install -c kotori_y scopy

Install with pypi
~~~~~~~~~~~~~~~~~~
>>> pip install cbdd-scopy

Molecular Pretreater
---------------------
The check and preparation for molecular structures are the necessary prerequisites for subsequent data analysis, especially for molecular resources downloaded from web sources. Considering its importance, the Scopy library provides `ScoPretreat`_ module to realize molecular preparation.

The `ScoPretreat`_ proivdes the following functions:

- Normalization of functional groups to a consistent format.
- Recombination of separated charges.
- Breaking of bonds to metal atoms.
- Competitive reionization to ensure strongest acids ionize first in partially ionize molecules.
- Tautomer enumeration and canonicalization.
- Neutralization of charges.
- Standardization or removal of stereochemistry information.
- Filtering of salt and solvent fragments.
- Generation of fragment, isotope, charge, tautomer or stereochemistry insensitive parent structures.
- Validations to identify molecules with unusual and potentially troublesome characteristics.

The functions can be customized by setting corresponding parameters according to the job demand, like diconnect metal ion.

>>> from rdkit import Chem
>>> from scopy.ScoPretreat import pretreat
>>>	#
>>> mol = Chem.MolFromSmiles('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
>>> sdm = pretreat.StandardizeMol()
>>> mol = sdm.disconnect_metals(mol) # diconnect metal ion.
>>> Chem.MolToSmiles(mol, isomericSmiles=True)
O=C([O-])c1ccc(C[S+2]([O-])[O-])cc1.[Na+]

Alternatively, users can achieve all above preparation steps by using the `StandardSmi` function.

>>> stdsmi = pretreat.StandardSmi('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
>>> stdsmi
O=C([O-])c1ccc(C[S](=O)=O)cc1

Drug-likeness Filter
---------------------
Drug-likeness is a conception that rationalizes the influence of simple physicochemical properties to in vivo molecular behaviors, with particular respect to solubility, absorption, permeability, metabolic stability and transporting effects. The application of drug-likeness rules to database construction will help senior executives more effectively.

The `ScoDruglikeness`_ module provides the calculation of physicochemical properties and the screening drug-likeness rules. This module can calculate 42 physicochemical properties (39 basic molecular properties and 3 comprehensive molecular evaluation scores), and implement 15 drug-likeness rules (11 drug-likeness rules, 2 macro-cycle molecule rules and 2 building block rules). More details see `overview`_.

Calculating Physicochemical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The class `ScoDruglikeness.molproperty_Lib.PC_properties` provides the calculation of 42 physicochemical properties, including 39 basic molecular properties and 3 comprehensive molecular evaluation scores.

The calculation of the physicochemical properties of 50 molecules will be taken as an example.

>>> import os
>>> from rdkit import Chem
>>> from scopy.ScoConfig import DemoDir
>>> mols = Chem.SDMolSupplier(os.path.join(DemoDir, '50.sdf'))
>>> mols = [mol for mol in mols if mol]
>>> len(mols)
50

Users can calculate different properties separately.

>>> from scopy.ScoDruglikeness import PC_properties
>>>	
>>> props = PC_properties(mols, n_jobs=4) #4 processors used to do the computation
>>> MW = props.CalculateMolWeight() #Calculate molecular weight.
>>> MW[:5]
[256.07, 288.06, 182.08, 578.14, 592.16]
>>> QEDnone = props.CalculateQEDnone() #Calculate QED using unit weights.
>>> QEDnone[:5]
[0.67, 0.42, 0.26, 0.1, 0.11]
>>> SAscore = props.CalculateSAscore() #Calculate Synthetic Accessibility Score
>>> SAscore[:5]
[4.08, 4.49, 3.56, 4.52, 4.55]
>>> NPscore = props.CalculateNPscore() #Calculate Natural Product-likeness Score
>>> NPscore[:5]
[0.64, 0.72, 1.14, 1.93, 2.04]

Alternatively, user can calculate multiple properties simultaneously through `GetProperties` method.

>>> mu_props = props.GetProperties(items=['MW','Vol','SAscore']) #The molecular weight, volume and SAscore to be calulated
>>> type(mu_props)
dict
>>> mu_props['MW'][:5]
[256.07, 288.06, 182.08, 578.14, 592
>>> mu_props['Vol'][:5]
[259.03, 276.61, 165.07, 549.94, 567.24]

Scopy propvide funtion to calculate physicochemical properties of single molecule in :mod:`ScoDruglikeness.molproperty`.

Screening under Drug-likeness Rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The class :mod:`ScoDruglikeness.rulesfilter_Lib.PC_rules` provides the screening of drug-likeness rules. In current version, the module can implement 15 drug-likeness rules, including 11 drug-likeness rules, 2 macro-cycle molecule rules and 2 building block rules.

>>> from scopy.ScoDruglikeness import PC_rules
>>> 
>>> rules = PC_rules(mols, n_jobs=4)
>>> ro5 = rules.CheckLipinskiRule() #Check the molecule whether math the requirements of Lipinski's Rule.
>>> ro5[:5]
[{'Disposed': 'Accepted', 'nViolate': 0},
 {'Disposed': 'Accepted', 'nViolate': 0},
 {'Disposed': 'Accepted', 'nViolate': 1},
 {'Disposed': 'Rejected', 'nViolate': 3}, #The disposed is Rejetced since violate 3 limitations of Lipinski's Rule
 {'Disposed': 'Rejected', 'nViolate': 3}]

User can also obtain more detailed information about the screening result.

>>> rules = PC_rules(mols, n_jobs=4, detail=True)
>>> ro5 = rules.CheckLipinskiRule() #Check the molecule whether math the requirements of Lipinski's Rule.
>>> ro5[2:5]
[{'MW': 182.08,
  'logP': -3.59,
  'nHD': 6,
  'nHA': 6,
  'Disposed': 'Accepted',
  'nViolate': 1},
 {'MW': 578.14,
  'logP': 3.0,
  'nHD': 10,
  'nHA': 12,
  'Disposed': 'Rejected',
  'nViolate': 3},
 {'MW': 592.16,
  'logP': 3.3,
  'nHD': 9,
  'nHA': 12,
  'Disposed': 'Rejected',
  'nViolate': 3}]

Considering the expert experience and different requirements in practical applications, users can customize their own screening rules through `Check_CustomizeRule` function.

>>> prop_kws = {'MW':[None,500], 'logP':[None, 5], 'nHD':[None,5], 'nHA':[None,10], 'TPSA':[None,140]} #The customized rule: MW<=500, logP<=5, nHD<=5, nHA<=10, TPSA<=140
>>> 
>>> custom = rules.CheckCustomizeRule(prop_kws)
>>> custom[:3]
[{'MW': 256.07,
  'logP': 2.83,
  'nHD': 3,
  'nHA': 3,
  'TPSA': 73.49,
  'nViolate': 0,
  'VioProp': []},
 {'MW': 288.06,
  'logP': 2.79,
  'nHD': 5,
  'nHA': 5,
  'TPSA': 113.95,
  'nViolate': 0,
  'VioProp': []},
 {'MW': 182.08,
  'logP': -3.59,
  'nHD': 6,
  'nHA': 6,
  'TPSA': 121.38,
  'nViolate': 1,
  'VioProp': ['nHD']}]

Scopy provides the visualization function to position the value of the queried compound within the selected drug-likeness rule ranges, which provide a benchmark for molecular assessment. See: `ScoVisualize.pc_depict.RuleRadar`_ function.

.. figure:: /image/user_guide/mol_basci_rule.png
	:width: 400px
	:align: center

Scopy also propvide funtion to screening rules properties of single molecule in :mod:`ScoDruglikeness.rulesfilter`.

Frequent Hitter Filter
------------------------
Frequent hitters refer to compounds which are repetitively identified as active hits in many different and independent biological assays covering a wide range of targets. Frequent hitters can be roughly divided into two categories: (1) compounds that interfere with elements of the assay formats or techniques thus causing undesirable false positive results; and (2) promiscuous compounds that can bind to different target thus triggering adverse reactions and other safety issues.

The `ScoFH`_ module provides 8 substructure filters for screening different types of FHs, including 4 assay interference substructure filters and 4 promiscuous compound substructure filters. More Details see `overview`_.

Assay Interference Substructure Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Assay interferences refer to compounds that interfere with elements of the assay formats or techniques thus causing undesirable false positive results. Such compounds will seriously interfere with the progress of drug research. class :mod:`ScoFH.fh_filter.FHfilter` provides 4 assay interference substructure filters (AlphaScreen_FHs, Luciferase_Inhibitory, Chelating and Alarm_NMR Filter) for the screening of AlphaScreen detection interferences, spectroscopic interferences, chelators and chemical reactive compounds, respectively.

>>> from scopy.ScoFH import FHfilter
>>>
>>> Filter = FHfilter(mols, n_jobs=4)
>>> res = Filter.Check_Alarm_NMR() #Here, Alarm_NMR Filter be used for screening the molecule.
>>> res[:3]
[{'Disposed': 'Accepted', 'Endpoint': 'Alarm_NMR'},
 {'Disposed': 'Rejected', 'Endpoint': 'Alarm_NMR'}, #Tthe status is 'Rejected' meant failed the ALARM NMR rule, 
 {'Disposed': 'Accepted', 'Endpoint': 'Alarm_NMR'}]

User can also obtain more detailed information about screening result.

>>> Filter = FHfilter(mols, n_jobs=4, detail=True, showSMILES=True)
>>> res = Filter.Check_Alarm_NMR()
[{'SMILES': 'OC1=Cc2c(O)cc(O)cc2[OH+][C-]1c1ccccc1',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'Alarm_NMR'},
 {'SMILES': 'OC1=C[C-]2[OH+]C(c3ccc(O)c(O)c3)=C(O)C=C2C(O)=C1',
  'Disposed': 'Rejected',
  'MatchedAtoms': [((19, 17, 18, 13, 14, 15, 16, 20),),
   ((19, 17, 16, 15, 14, 13, 18), (20, 16, 15, 14, 13, 18, 17))],
  'MatchedNames': ['[OH]c1ccccc1O', 'c1ccccc1O'],
  'Endpoint': 'Alarm_NMR'},
 {'SMILES': 'OCC(O)C(O)C(O)C(O)CO',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'Alarm_NMR'}]

Promiscuous Compound Substructure Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The promiscuity is defined as the ability to specifically bind to different macro-molecular targets. These multiple interactions can include unintended targets, thus triggering adverse reactions and other safety issues. class :mod:`ScoFH.fh_filter.FHfilter` provides 4 frequently-used promiscuous compound substructure filters, such as PAINS, BMS Filter, AlphaScreen_GST_FHs and AlphaScreen_HIS_FHs.

>>> res = Filter.Check_PAINS(mol, detail=True) #Here, PAINS Filter used for screening the molecule.
>>> res[5:8]
[{'SMILES': 'O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
  'Disposed': 'Rejected',
  'MatchedAtoms': [((13, 14, 15, 16, 17, 18, 19, 20),)],
  'MatchedNames': ['Catechol_A'],
  'Endpoint': 'Pains'},
 {'SMILES': 'O=c1cc(-c2ccc(O)cc2)oc2cc(OC3OC(CO)C(O)C(O)C3O)cc(O)c12',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'Pains'},
 {'SMILES': 'O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'Pains'}]

`Scovisualize.highlight.HighlightAtoms` allows user to conduct further analysis and molecular optimization, which also provides intuitive information about the vigilant alerts.

.. figure:: /image/user_guide/PAINS.svg
	:width: 400px
	:align: center

Toxicity Filter
----------------
Toxicity refers to the measure of poisonous or toxic effect on an organ or a whole organism. Toxicity is one of the main reasons for attrition in the drug development process. It is reported that more than 15% of new approved FDA chemical entitles (between 1975 and 2009) have received more than once black-box warnings, and some of them have been withdrawn from the market due to the toxicity and safety issues. In addition, the requirements for molecular safety are not only limited to the human beings. The environmental influence of drugs has also aroused great concern.

`ScoTox`_ module provides 11 toxicophore filters, including 5 human related toxicity substructure filters, 3 environment related toxicity substructure filters and 3 comprehensive substructure filters. More details see: `overview`_.

Human Toxic Compound Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~
For toxicity to human beings, 5 credible toxicophore filters are used to evaluate the potential toxicity of query compounds, from broad toxicity and acute toxicity, to carcinogenicity and mutagenicity.

>>> from scopy.ScoTox import Toxfilter
>>> Filter = Toxfilter(mols, detail=True, showSMILES=True)
>>> res = Filter.Check_Genotoxic_Carcinogenicity_Mutagenicity() #This Filter related with carcinogenicity and mutagenicity.
>>> res[:3]
[{'SMILES': 'OC1=Cc2c(O)cc(O)cc2[OH+][C-]1c1ccccc1',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'Genotoxic_Carcinogenicity_Mutagenicity'},
 {'SMILES': 'OC1=C[C-]2[OH+]C(c3ccc(O)c(O)c3)=C(O)C=C2C(O)=C1',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'Genotoxic_Carcinogenicity_Mutagenicity'},
 {'SMILES': 'OCC(O)C(O)C(O)C(O)CO',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'Genotoxic_Carcinogenicity_Mutagenicity'}]

Environmental Toxic Compound Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Given the stringent requirements for environmental safety, the Scopy library provides 3 substructure filters for the evaluation of molecular biodegradability and potential aquatic toxicity.

>>> res = Filter.Check_NonBiodegradable()
>>> res[:3]
[{'SMILES': 'OC1=Cc2c(O)cc(O)cc2[OH+][C-]1c1ccccc1',
  'Disposed': 'Rejected',
  'MatchedAtoms': [((1,), (15,))],
  'MatchedNames': ['MoreThanTwoHydroxyOnAromaticRing'],
  'Endpoint': 'NonBiodegradable'},
 {'SMILES': 'OC1=C[C-]2[OH+]C(c3ccc(O)c(O)c3)=C(O)C=C2C(O)=C1',
  'Disposed': 'Rejected',
  'MatchedAtoms': [((16,), (17,))],
  'MatchedNames': ['MoreThanTwoHydroxyOnAromaticRing'],
  'Endpoint': 'NonBiodegradable'},
 {'SMILES': 'OCC(O)C(O)C(O)C(O)CO',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'NonBiodegradable'}]

Comprehensive Toxic Compound Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To simplify screening process and draw lessons from existing screening tools, the Scopy library has integrated 3 comprehensive filters from FAF-Drugs4, SureChEMBL and Brenk et.al work.

>>> res = Filter.Check_SureChEMBL()
>>> res[:3]
[{'SMILES': 'OC1=Cc2c(O)cc(O)cc2[OH+][C-]1c1ccccc1',
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'SureChEMBL'},
 {'SMILES': 'OC1=C[C-]2[OH+]C(c3ccc(O)c(O)c3)=C(O)C=C2C(O)=C1',
  'Disposed': 'Rejected',
  'MatchedAtoms': [((2, 1, 10, 9, 8, 7, 6, 5),)],
  'MatchedNames': ['polyenes'],
  'Endpoint': 'SureChEMBL'},
 {'SMILES': 'OCC(O)C(O)C(O)C(O)CO',      
  'Disposed': 'Accepted',
  'MatchedAtoms': ['-'],
  'MatchedNames': ['-'],
  'Endpoint': 'SureChEMBL'}]

Chemical Space Exploer
------------------------
A desirable database is demanded to own wide chemical space, which will greatly benefits the efficiency and success rate of drug development. To analyze the chemical diversity of screening databases, the `ScoRepresent`_ can calculate 2 molecular scaffolds, 6 substructure descriptors and 2 fingerprints.

Framework Calculation
~~~~~~~~~~~~~~~~~~~~~~~
The functions from `ScoRepresent.scaffolds` can calculate molecular Murcko scaffold and carbon skeleton and summarize the number of scaffold occurrence in the database. Then the data can be used to generate the cloud gram via `ScoVisualize.mcloud.ShowMcloud`_ function. 

The function `ScoRepresent.scaffolds.CountMurckoFramework` can calculate the Murcko framework and count the frequency of corresponding frameworks.

>>> from scopy.ScoRepresent import CountMurckoFramework
>>>		
>>> scount = CountScaffold(mols)
>>> type(scount)
>>> dict
>>> len(scount)
>>> 50
>>> list(scount.keys())[:3]
['C1=C[C-](c2ccccc2)[OH+]c2ccccc21',
 'C1=CC2=CC=C(c3ccccc3)[OH+][C-]2C=C1',
 'c1ccc(C2CC(c3cccc4c3OC(c3ccccc3)CC4)c3ccccc3O2)cc1']
>>> list(scount.values())[:3]
[1, 1, 2]

.. figure:: /image/user_guide/mcloud.png
	:width: 500px
	:align: center
	

Fingerprint Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~
Different type of substructure descriptors can characterize the structural features of the studied molecules from different viewpoints.

The `ScoRepresent`_ module provides the calculation of 6 descriptors (MACCS, EFG, IFG, EState, GhoseCrippen and PubChem) and 2 fingerprints (Morgan Family and Daylight Fingerprint). More Details see `overview`_. As for Morgan Family, 2 and 1024 chosen as the default radius and the number of bits. Besides, minimum and maximum distance for Daylight fingerprint, whose default the number of bits is 2048, set default as 1 and 7.

>>> from scopy.ScoRepresent import CalculateEFG
>>> fps = CalculateEFG(mols, useCount=True, n_jobs=4)
>>> type(fps)
numpy.ndarray
>>> fps.shape
(50, 583)

Screening Visualizer
--------------------
In the case of early drug discovery, data visualized as a gram or diagram can provide a simplified view of multidimensional property and ideally reveal correlations. The `ScoVisualize`_ module provides four different visualization functions, including basic feature radar charts, feature-feature related scatter diagram, functional group marker gram and cloud gram.

PC Visualizer
~~~~~~~~~~~~~
The submodule `ScoVisualize.pc_depict` module provides the visualization of PC properties distribution and drug-likeness rules.

Proprty Matrix
""""""""""""""
The proprty matrix (feature-feature related scatter diagram) can present the correlation between different features and assessment score.

>>> from scopy.ScoVisualize import pc_depict
>>> 
>>> fig = pc_depict.PropMatrix(mols)
>>> fig
<Figure size 1567x989 with 36 Axes>

.. figure:: /image/user_guide/50_matrix.png
	:width: 600px
	:align: center

	The matrix of logP, TPSA, MW, nRot, nHD and nHA

Default properties of matrix are logP, TPSA, MW, nRot, nHD and nHA. Users can customize their own features.

>>> fig = pc_depict.PropMatrix(mols, n_jobs=4, items=['MW', 'Vol', 'Dense']) #Mw, Vol and Dense to be shown.

.. figure:: /image/user_guide/50_matrix_2.png
	:width: 500px
	:align: center

	The matrix of MW, Vol and Dense

Basic Property Radar
""""""""""""""""""""
The radar chart can be used to visualize the physicochemical properties for a single compound or a whole database. In addition, user can position the value of the queried compound within the selected drug-likeness rule range, which may provide a benchmark for molecular assessment.

>>> fig = pc_depict.RuleRadar(mols[0])
>>> fig
<Figure size 640x480 with 1 Axes>

.. figure:: /image/user_guide/mol_basci_rule.png
	:width: 500px
	:align: center

Fragment visualizer
~~~~~~~~~~~~~~~~~~~
The function `Scovisualize.highlight.HighlightAtoms` can highlight the flagged substructures, which help user to conduct further analysis.

>>> from scopy.ScoVisualize import highlight
>>> 
>>> fig = highlight.HighlightAtoms(mols[0], highlightAtoms=[13, 14, 15, 16, 17, 18, 19, 20]) #highlightAtoms obtained from function Check_PAINS()
>>> type(fig)
IPython.core.display.SVG

.. figure:: /image/user_guide/PAINS.svg
	:width: 500px
	:align: center

Framework Visualizer
~~~~~~~~~~~~~~~~~~~~~~
The function `ScoVisualize.mcloud.ShowMcloud` can help the evaluation of database diversity and structure characteristics. The frequency of specific scaffold is indicated by the size of the respective structural image. With the application of cloud gram, users can easily explore the top-ranked scaffolds and the whole chemical space of the screening database.

.. note::
	This module should run under a Java environment and the script retrived from `Peter Ertl`_

>>> scaffolds = os.path.join(ScoConfig.DemoDir, 'scaffolds.txt') #The file storing the frameworks and corresponding frequency.
>>> mcloud.ShowMcloud(file=scaffolds, number=200, skip=1) #The skip parameter is used to skip the most frequent framework (here skipping benzene ring).

.. figure:: /image/user_guide/mcloud.png
	:width: 500px
	:align: center
	

.. _`ScoPretreat`: ./modules/scopy.ScoPretreat.html
.. _`ScoDruglikeness`: ./modules/scopy.ScoDruglikeness.html
.. _`ScoFH`: ./modules/scopy.ScoFH.html
.. _`ScoTox`: ./modules/scopy.ScoTox.html
.. _`ScoRepresent`: ./modules/scopy.ScoRepresent.html
.. _`ScoVisualize`: ./modules/scopy.ScoVisualize.html
.. _`overview`: ./overview.html#feature-overview
.. _`Scovisualize.highlight.HighlightAtoms`: #fragment-visualizer
.. _`ScoVisualize.pc_depict.	RuleRadar`: #basic-property-radar
.. _`ScoVisualize.mcloud.ShowMcloud`: #framework-visualizer
.. _`Frequent Hitters Filter`: #frequent-hitter-filter
.. _`Toxicity Filter`: #toxicity-filter
.. _`Peter Ertl`: https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-12