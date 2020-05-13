..  -*- coding: utf-8 -*-

Getting Started with Scopy
==========================
This document intends to provide users with the basic operation methods of Scopy. If you find any mistake or have suggestions for improvements, please either fix them in the source document (the .py file) or send to the mailing list: oriental-cds@163.com and kotori@cbdd.me.

Installing the Scopy package
-----------------------------
Scopy has been successfully tested on Linux and Windows systems under python3 enviroment.

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
The check and preparation for molecular structures is the necessary prerequisite for subsequent data analysis, especially for molecular resources downloaded from web sources. Considering its importance, the Scopy library provides scopy.pretreat.pretreat module to realize molecular preparation.

The :mod:`scopy.pretreat.pretreat` proivdes the following functions:

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

Users can diconnect metal ion.

>>> from rdkit import Chem
>>> from scopy.pretreat import pretreat
>>>	
>>> mol = Chem.MolFromSmiles('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
>>> sdm = pretreat.StandardizeMol()
>>> mol = sdm.disconnect_metals(mol) # diconnect metal ion.
>>> Chem.MolToSmiles(mol, isomericSmiles=True)
O=C([O-])c1ccc(C[S+2]([O-])[O-])cc1.[Na+]

And pretreat the molecular structure using all functions.

>>> stdsmi = pretreat.StandardSmi('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
>>> stdsmi
O=C([O-])c1ccc(C[S](=O)=O)cc1

Drug-likeness Filter
---------------------
Drug-likeness is a conception that rationalizes the influence of simple physicochemical properties to in vivo molecular behaviors, with particular respect to solubility, absorption, permeability, metabolic stability and transporting effects. The application of drug-likeness rules to database construction will help senior executives more effectively. The Drug-likeness Filter is implemented in the :mod:`scopy.druglikeness` package.

The :mod:`scopy.druglikeness` package provides the calculation of physicochemical properties and the screening drug-likeness rules. :mod:`scopy.druglikeness` package can calculate 42 physicochemical properties (39 basic molecular properties and 3 comprehensive molecular evaluation scores), and implement 15 drug-likeness rules (11 drug-likeness rules, 2 macro-cycle molecule rules and 2 building block rules). More details see `overview`_.

Calculating Physicochemical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`scopy.druglikeness.molproperty` module provides the calculation of 42 physicochemical properties, including 39 basic molecular properties and 3 comprehensive molecular evaluation scores.

>>> mol = Chem.MolFromSmiles('Cc1cc(O)cc(N=C2C=CC(=O)C(O)=C2)c1')
>>> mol
<rdkit.Chem.rdchem.Mol object at 0x0000020879B4E120>

.. figure:: /image/user_guide/demo_mol.svg
	:width: 300px
	:align: center

	The molecule used as the example in this document.

Users can calculate different properties separately.

>>> from scopy.druglikeness import molproperty
>>>	
>>> MW = molproperty.CalculateMolWeight(mol) #Calculate molecular weight.
>>> MW
229.07
>>> QEDnone = molproperty.CalculateQEDnone(mol) #Calculate QED using unit weights.
>>> QEDnone
0.79
>>> SAscore = molproperty.CalculateSAscore(mol) #Calculate Synthetic Accessibility Score
>>> SAscore
2.96
>>> NPscore = molproperty.CalculateNPscore(mol) #Calculate Natural Product-likeness Score
>>> NPscore
0.49

Besides, user can also calculate different property simultaneously through `molproperty.GetProperties` function.

>>> props = molproperty.GetProperties(mol, items=['MW','Vol','SAscore']) #The molecular weight, volume and SAscore to be calulated
>>> props
{'MW': 229.07, 'Vol': 235.2, 'SAscore': 2.96}

When user needs to calculte properties of a set of molecules, :mod:`scopy.druglikeness.druglikeness` module can be used for the fast implementaiton.

>>> import os
>>> from scopy import ScoConfig
>>> from scopy.druglikeness import druglikeness
>>> suppl = Chem.SDMolSupplier(os.path.join(ScoConfig.DemoDir, '760.sdf'))
>>> mols = [mol for mol in suppl]
>>>	
>>> pc = druglikeness.PC_properties(mols=mols, n_jobs=4) #4 processors used to do the computation.
>>> res = pc.CalculateMolWeight()
>>> len(res)
760
>>> type(res)
list
>>> res[:10]
[256.07, 288.06, 182.08, 578.14, 592.16, 286.05, 432.11, 270.05, 448.1, 578.16]

Screening under Drug-likeness Rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`scopy.druglikeness.rulesfilter` module provides the screening of drug-likeness rules. In current version, the module can implement 15 drug-likeness rules, including 11 drug-likeness rules, 2 macro-cycle molecule rules and 2 building block rules.

>>> from scopy.druglikeness import rulesfilter
>>> res = rulesfilter.CheckLipinskiRule(mol) #Check the molecule whether math the requirements of Lipinski's Rule.
>>> res
{'Disposed': 'Accepted', 'nViolate': 0}

In above example, the molecule do does not violate any property limited limitations in of Lipinski's Rule. Thus its status is 'Accepted'.

Besides, users can obtain more detailed information about the screening result.

>>> res = rulesfilter.CheckLipinskiRule(mol, detail=True, showSMILES=True)
>>> res
{'SMILES': 'Cc1cc(O)cc(N=C2C=CC(=O)C(O)=C2)c1',
 'MW': 229.07,
 'logP': 2.35,
 'nHD': 2,
 'nHA': 4,
 'Disposed': 'Accepted',
 'nViolate': 0}

Considering the expert experience and different requirements in practical applications, users can customize their own screening rules through `rulesfilter.Check_CustomizeRule` function.

>>> prop_kws = {'MW':[None,500], 'logP':[None, 5], 'nHD':[None,5], 
... 			'nHA':[None,10], 'TPSA':[None,140]} #The customized rule: MW<=500, logP<=5, nHD<=5, nHA<=10, TPSA<=140
>>> res = rulesfilter.Check_CustomizeRule(mol, prop_kws=prop_kws, detail=True)
>>> res
{'MW': 229.07,
 'logP': 2.35,
 'nHD': 2,
 'nHA': 4,
 'TPSA': 69.89,
 'nViolate': 0,
 'VioProp': []}

Scopy provides the visualization function to position the value of the queried compound within the selected drug-likeness rule ranges, which provide a benchmark for molecular assessment. See: `visualize.rule_radar`_ function.

.. figure:: /image/user_guide/radar_1.png
	:width: 400px
	:align: center

Similarly, :mod:`scopy.druglikeness.druglikeness` module can be used to evaluate the potential of a group of molecules.

>>> rule = druglikeness.PC_rules(mols, detail=True, n_jobs=4) #4 processors used to do the computation.
>>> res = rule.CheckLipinskiRule()
>>> len(res)
760
>>> type(res)
list
>>> res[:3]
[{'MW': 256.26, 'logP': 2.83, 'nHD': 3, 'nHA': 3, 'Disposed': 'Accepted', 'nViolate': 0},
 {'MW': 288.25, 'logP': 2.79, 'nHD': 5, 'nHA': 5, 'Disposed': 'Accepted', 'nViolate': 0},
 {'MW': 182.17, 'logP': -3.59, 'nHD': 6, 'nHA': 6, 'Disposed': 'Accepted', 'nViolate': 1}]

Frequent hitter Filter
------------------------
Frequent hitters refer to compounds which are repetitively identified as active hits in many different and independent biological assays covering a wide range of targets. Frequent hitters can be roughly divided into two categories: (1) compounds that interfere with elements of the assay formats or techniques thus causing undesirable false positive results; and (2) promiscuous compounds that can bind to different target thus triggering adverse reactions and other safety issues.

The :mod:`scopy.structure_alert` package provides 8 substructure filters for screening different types of FHs, including 4 assay interference substructure filters and 4 promiscuous compound substructure filters. More Details see `overview`_.

Assay Interference Substructure Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Assay interferences refer to compounds that interfere with elements of the assay formats or techniques thus causing undesirable false positive results. Such compounds will seriously interfere with the progress of drug research. :mod:`scopy.structure_alert.FilterWithSmarts` module provides 4 assay interference substructure filters (AlphaScreen_FHs, Luciferase_Inhibitory, Chelating and Alarm_NMR Filter) for the screening of AlphaScreen detection interferences, spectroscopic interferences, chelators and chemical reactive compounds, respectively.

>>> from scopy.structure_alert import FilterWithSmarts
>>> res = FilterWithSmarts.Check_Alarm_NMR(mol) #Here, Alarm_NMR Filter be used for screening the molecule.
>>> res
{'Disposed': 'Rejected', 'Endpoint': 'Alarm_NMR'}

In the above example, the molecule failed the ALARM NMR rule.

User can also obtain more detailed information about screening result.

>>> res = FilterWithSmarts.Check_Alarm_NMR(mol, detail=True, showSMILES=True)
>>> res
{'SMILES': 'Cc1cc(O)cc(N=C2C=CC(=O)C(O)=C2)c1',
 'Disposed': 'Rejected',
 'MatchedAtoms': [((9, 10, 11, 12, 13), (15, 13, 11, 12, 10)),
  ((4, 3, 5, 6, 7, 16, 1, 2),),
  ((7, 6, 5, 3, 2, 1, 16),),
  ((4, 3, 2, 1, 16, 6, 5),)],
 'MatchedNames': ['C=CC(=O)C', '[OH]c1cc(N)ccc1', 'c1ccccc1N', 'c1ccccc1O'],
 'Endpoint': 'Alarm_NMR'}

Promiscuous Compound Substructure Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The promiscuity is defined as the ability to specifically bind to different macro-molecular targets. These multiple interactions can include unintended targets, thus triggering adverse reactions and other safety issues. :mod:`scopy.structure_alert.FilterWithSmarts` module provides 4 frequently-used promiscuous compound substructure filters, such as PAINS, BMS Filter, AlphaScreen_GST_FHs and AlphaScreen_HIS_FHs.

>>> res = FilterWithSmarts.Check_PAINS(mol, detail=True) #Here, PAINS Filter used for screening the molecule.
>>> res
{'Disposed': 'Rejected',
 'MatchedAtoms': [((7, 8, 9, 10, 11, 12, 13, 15),)],
 'MatchedNames': ['Quinone_A'],
 'Endpoint': 'Pains'}

By applying `visualize.HighlightAtoms.highlight`_ function, user conduct further analysis and molecular optimization, which also provide intuitive information about the vigilant alerts.

.. figure:: /image/user_guide/highlight_1.svg
	:width: 400px
	:align: center

Toxicity Filter
----------------
Toxicity refers to the measure of poisonous or toxic effect on an organ or a whole organism. Toxicity is one of the main reasons for attrition in the drug development process. It is reported that more than 15% of new approved FDA chemical entitles (between 1975 and 2009) have received more than once black-box warnings, and some of them have been withdrawn from the market due to the toxicity and safety issues. In addition, the requirements for molecular safety are not only limited to the human beings. The environmental influence of drugs has also aroused great concern. 

:mod:`scopy.structure_alert` package provides 11 toxicophore filters, including 5 human related toxicity substructure filters, 3 environment related toxicity substructure filters and 3 comprehensive substructure filters. More details see: `overview`_.

Human Toxic Compound Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~
For toxicity to human beings, 5 credible toxicophore filters are used to evaluate the potential toxicity of query compounds, from broad toxicity and acute toxicity, to carcinogenicity and mutagenicity.

>>> res = FilterWithSmarts.Check_Genotoxic_Carcinogenicity_Mutagenicity(mol) #This Filter related with carcinogenicity and mutagenicity.
>>> res
{'Disposed': 'Rejected', 'Endpoint': 'Genotoxic_Carcinogenicity_Mutagenicity'}

User can also check more detailed information.

>>> res = FilterWithSmarts.Check_Genotoxic_Carcinogenicity_Mutagenicity(mol, detail=True)
>>> res
{'Disposed': 'Rejected',
 'MatchedAtoms': [((1, 0, 15, 16, 13), (12, 13, 15, 16, 0))], #It means there two corresponding substructure in this molecule
 'MatchedNames': ['α, β-Unsaturated carbonyls'],
 'Endpoint': 'Genotoxic_Carcinogenicity_Mutagenicity'}

Environmental Toxic Compound Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Given the stringent requirements for environmental safety, the Scopy library provides 3 substructure filters for the evaluation of molecular biodegradability and potential aquatic toxicity.

>>> res = FilterWithSmarts.Check_NonBiodegradable(pains_mol, detail=True)
>>> res
{'Disposed': 'Rejected', 'Endpoint': 'NonBiodegradable'}

Comprehensive Toxic Compound Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To simplify screening process and draw lessons from existing screening tools, the Scopy library has integrated 3 comprehensive filters from FAF-Drugs4, SureChEMBL and Brenk et.al work.

>>> res = FilterWithSmarts.Check_SureChEMBL(mol)
{'Disposed': 'Accepted', 'Endpoint': 'SureChEMBL'}

Multiprocessing Filter
-----------------------
The :mod:`scopy.structure_alert.SmartsFilter` module provides the tool to screen molecule library under `Frequent Hitters Filter`_ and (or) `Toxicity Filter`_ with multiprocess technology.

>>> from scopy.structure_alert import SmartsFilter
>>>		
>>> Screener = SmartsFilter.Filter(mols, detail=True, n_jobs=4) #4 processors used to do the screening.
>>> res = Screener.Check_PAINS() #Here, PAINS Filter used for screening library.
>>> type(res)
list
>>> len(res)
760
>>>	
>>> res[0]
{'Disposed': 'Accepted', 'MatchedAtoms': ['-'], 'MatchedNames': ['-'], 'Endpoint': 'Pains'}
>>> res[207]
{'Disposed': 'Rejected', 'MatchedAtoms': [((7, 16, 15, 17, 18, 19, 20, 21, 14),)], 'MatchedNames': ['Mannich_A'], 'Endpoint': 'Pains'}

Chemical Space Exploer
------------------------
A desirable database is demanded to own wide chemical space, which will greatly benefits the efficiency and success rate of drug development. To analyze the chemical diversity of screening databases, the Scopy library designs a special module for the calculation of 2 molecular scaffolds, 6 substructure descriptors and 2 fingerprints. 

Framework Calculation
~~~~~~~~~~~~~~~~~~~~~~~
The function `mcloud.CountScaffold` can calculate molecular Murcko scaffold and carbon skeleton and summarize the number of scaffold occurrence in the database. Then the data can be used to generate the cloud gram via `visualize.mcloud`_ function. 

The function `mcloud.CountScaffold` can calculate the framework and count the frequency of corresponding frameworks.

>>> from scopy.visualize import mcloud
>>>		
>>> scount = mcloud.CountScaffold(mols)
>>> type(scount)
>>> dict
>>> len(scount)
>>> 760
>>> list(scount.keys())[:3]
['C1=C[C-](c2ccccc2)[OH+]c2ccccc21',
 'C1=CC2=CC=C(c3ccccc3)[OH+][C-]2C=C1',
 'c1ccc(C2CC(c3cccc4c3OC(c3ccccc3)CC4)c3ccccc3O2)cc1']
>>> list(scount.values())[:3]
[1, 3, 3]

.. figure:: /image/user_guide/mcloud_1.png
	:width: 500px
	:align: center
	

Fingerprint Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~
With different definitions, fingerprints (descriptors) can characterize molecules from different angles. Through calculating similarity or distance among molecular fingerprints (descriptors), the spatial density of compound libraries can be evaluated.

The :mod:`scopy.fingerprint` package provides the calculation of 6 descriptors (MACCS, EFG, IFG, EState, GhoseCrippen and PubChem) and 2 fingerprints (Morgan Family and Daylight Fingerprint). More Details see `overview`_. As for Morgan Family, 2 and 1024 chosen as the default radius and the number of bits. Besides, minimum and maximum distance for Daylight fingerprint, whose default the number of bits is 2048, set default as 1 and 7.

>>> from scopy.fingerprint import fingerprints
>>> fps = fingerprints.CalculateEFG(mols, useCount=True, n_jobs=4)
>>> fps.shape
(760, 583)

Screening Visualizer
--------------------
In the case of early drug discovery, data visualized as a gram or diagram can provide a simplified view of multidimensional property and ideally reveal correlations. The :mod:`scopy.visualize` module provides four different visualization functions, including basic feature radar charts, feature-feature related scatter diagram, functional group marker gram and cloud gram.

PC Visualizer
~~~~~~~~~~~~~
The :mod:`scopy.visualize.pc_depict` module provides the visualization of PC properties distribution and drug-likeness rules.

Proprty Matrix
""""""""""""""
The proprty matrix (feature-feature related scatter diagram) can present the correlation between different features and assessment score.

>>> from scopy.visualize import pc_depict
>>> #
>>> fig = pc_depict.prop_matrix()
>>> fig
<Figure size 1567x989 with 36 Axes>

.. figure:: /image/760_matrix.png
	:width: 500px
	:align: center

	The matrix of logP, TPSA, MW, nRot, nHD and nHA

Default properties of matrix are logP, TPSA, MW, nRot, nHD and nHA. Users can customize their own features.

>>> fig = pc_depict.prop_matrix(mols, n_jobs=4, items=['MW', 'Vol', 'Dense']) #Mw, Vol and Dense to be shown.

.. figure:: /image/760_matrix_2.png
	:width: 500px
	:align: center

	The matrix of MW, Vol and Dense

Basic Property Radar
""""""""""""""""""""
The radar chart can be used to position the value of the queried compound within the selected drug-likeness rule ranges, which provide a benchmark for molecular assessment.

>>> fig = pc_depict.rule_radar(mol)
>>> fig
<Figure size 640x480 with 1 Axes>

.. figure:: /image/mol_basci_rule.png
	:width: 500px
	:align: center

Fragment visualizer
~~~~~~~~~~~~~~~~~~~
The :mod:`scopy.visualize.highlight` module can highlight the flagged substructures, which help user to conduct further analysis.

>>> from scopy.visualize import highlight
>>> #
>>> fig = highlight.HighlightAtoms(mol, highlightAtoms=[7, 8, 9, 10, 11, 12, 13, 15]) #highlightAtoms obtained from function Check_PAINS()
>>> type(fig)
IPython.core.display.SVG

.. figure:: /image/user_guide/highlight_1.svg
	:width: 500px
	:align: center

Framework Visualizer
~~~~~~~~~~~~~~~~~~~~~~
The function `mcloud.ShowMcloud` can help the evaluation of database diversity and structure characteristics. The frequency of specific scaffold is indicated by the size of the respective structural image. With the application of cloud gram, users can easily explore the top-ranked scaffolds and the whole chemical space of the screening database.

.. note::
	This module should run under a Java environment and the script retrived from `Peter Ertl`_

>>> scaffolds = os.path.join(ScoConfig.DemoDir, 'scaffolds.txt') #The file storing the frameworks and corresponding frequency.
>>> mcloud.ShowMcloud(file=scaffolds, number=200, skip=1) #The skip parameter is used to skip the most frequent framework (here skipping benzene ring).

.. figure:: /image/user_guide/mcloud_1.png
	:width: 500px
	:align: center
	


.. _`overview`: ./overview.html#feature-overview
.. _`visualize.HighlightAtoms.highlight`: #fragment-visualizer
.. _`visualize.rule_radar`: #basic-property-radar
.. _`visualize.mcloud`: #framework-visualizer
.. _`Frequent Hitters Filter`: #frequent-hitter-filter
.. _`Toxicity Filter`: #toxicity-filter
.. _`Peter Ertl`: https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-12