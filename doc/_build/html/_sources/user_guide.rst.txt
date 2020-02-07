..  -*- coding: utf-8 -*-

Getting Started with Scopy
==========================
This document is intended to provide an overview of how one can use the Scopy functionality from Python. If you find mistakes, or have suggestions for improvements, please either fix them yourselves in the source document (the .py file) or send them to the mailing list: oriental-cds@163.com and yzjkid9@gmail.com.

Installing the Scopy package
-----------------------------
PyBioMed has been successfully tested on Linux and Windows systems under python3 enviroment.

Dependencies
~~~~~~~~~~~~
.. code-block:: python

	RDkit>=2019.03.1
	Numpy
	Matplotlib

Install from source
~~~~~~~~~~~~~~~~~~~
>>> git clone git@github.com:kotori-y/Scopy.git && cd scopy
>>> [sudo] python setup.py install

Molecular Pretreater
---------------------
The :mod:`scopy.pretreat` can pretreat the molecular structure. The :mod:`scopy.pretreat` proivdes the following functions:

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

>>> from scopy.pretreat import pretreat
>>> mol = Chem.MolFromSmiles('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
>>> sdm = StandardizeMol()
>>> mol = sdm.disconnect_metals(mol)
>>> Chem.MolToSmiles(mol, isomericSmiles=True)
O=C([O-])c1ccc(C[S+2]([O-])[O-])cc1.[Na+]

Drug-likeness Filter
--------------------
The :mod:`druglikeness` provides the method to analyse the physicochemical (PC) properties and filter compounds based on PC-derived rules. 

Calculating PC properties
~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`druglikeness.molproperty` module provides the tool to calculate PC properties.

>>> import os
>>> from rdkit import Chem
>>> from scopy import ScoConfig
>>> from scopy.druglikeness import molproperty

>>> mol = Chem.MolFromMolFile(os.path.join(ScoConfig.DemoDir,'mol.sdf'))
>>> mol
<rdkit.Chem.rdchem.Mol object at 0x0000020879B4E120>

User could calculate different properties respectively

>>> MW = molproperty.CalculateMolWeight(mol)
>>> MW
229.24
>>> logP = molproperty.CalculateLogP(mol)
>>> logP
2.35
>>> nHD = molproperty.CalculateNumHDonors(mol)
>>> nHD
2

Beside these basic properties, some spcecial properties could also be calculated

>>> QEDnone = molproperty.CalculateQEDnone(mol)
>>> QEDnone
0.79

QED (quantitative estimate of drug-likeness) is a measure of drug-likeness. More datails: `Nat Chem 2012`_

>>> SAscore = molproperty.CalculateSAscore(mol)
>>> SAscore
2.96

SA (Synthetic Accessibility) score measure the synthetic accessibility of a molecule based on molecular complexity and fragment contributions. More details: `J Cheminform 2009`_

>>> NPscore = molproperty.CalculateNPscore(mol)
>>> NPscore
0.49

NP (Natural Product-likeness) score measure the natural product-likeness of a molecule. More details: `J Chem Inf Model 2008`_

.. _`Nat Chem 2012`: https://www.nature.com/nchem/journal/v4/n2/abs/nchem.1243.html
.. _`J Cheminform 2009`: https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8
.. _`J Chem Inf Model 2008`: https://pubs.acs.org/doi/abs/10.1021/ci700286x

User could also calculate multi-property at once through :mod:`molproperty.GetProperties`.

>>> props = molproperty.GetProperties(mol, items=['MW','Vol','SAscore'])
>>> props
{'MW': 229.07, 'Vol': 235.2, 'SAscore': 2.96}

The function return a `dict`, user could pass properties need to be calculated to parameter `item`, defaults to the whole (45) properties.

When calculating the property of multiple molecules, in addition to repeatedly calling the function in :mod:`druglikeness.molproperty`, you can also use :mod:`druglikeness.druglikeness` module, which is more time-saveing since using multiprocessing.

>>> from scopy.druglikeness import druglikeness
>>> suppl = Chem.SDMolSupplier(os.path.join(ScoConfig.DemoDir, '760.sdf'))
>>> mols = (mol for mol in suppl if mol)

>>> pc = druglikeness.PC_properties(mols=mols, n_jobs=4)
>>> res = pc.CalculateMolWeight()
>>> len(res)
760
>>> type(res)
<class 'list'>
>>> res[:10]
[256.26, 288.25, 182.17, 578.53, 592.55, 286.24, 432.38, 270.24, 448.38, 578.52]

The function return a `list`. Parameter `mols` should be an iterable object (i.g. `list`, `tuple` or `generator`) and `n_jobs` is the number of CPUs to use to do the computation, -1 means using all processors.

Filtering molecule under PC-derived rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`druglikeness.rulesfilter` module provides the tool to analyse PC properties

>>> from scopy.druglikeness import rulesfilter
>>> res = rulesfilter.CheckLipinskiRule(mol)
>>> res
{'Disposed': 'Accepted', 'nViolate': 0}

The function return a `dict`, the field :mod:`Disposed` represents compound state after filter applied (**Rejected** meant the compound rejected by filter, **Accepted** for accepted); :mod:`nViolate` represents the number of PC property violated by compound.

In above example, the compound do not violate any property limited in Lipinski Rule thus its status is 'Accepted'.

Besides, the specific value of each propety would be returned if the :mod:`detail` has been set as :mod:`True` and the SMILES would be also returned if the :mod:`showSMILES` has been set as :mod:`True`.

>>> res = rulesfilter.CheckLipinskiRule(mol, detail=True, showSMILES=True)
>>> res
{'MW': 229.24, 'logP': 2.35, 'nHD': 2, 'nHA': 4, 'Disposed': 'Accepted', 'nViolate': 0, 'SMILES': 'Cc1cc(O)cc(/N=C2/C=CC(=O)C(O)=C2)c1'}

You also could customize the filter by your experience

>>> prop_kws = {'MW':[100,500], 'nHB':[5,10], 'QEDmean':[0.8,None]}
>>> res = rulesfilter.Check_CustomizeRule(mol, prop_kws=prop_kws, detail=True)
>>> res
{'MW': 229.24, 'nHB': 6, 'QEDmean': 0.73, 'nViolate': 1, 'VioProp': ['QEDmean']}

The customize rule should be a `dict`, key of `dict` is abbreviation name of property and value is the limited range.

Samely, :mod:`druglikeness.druglikeness` could also be used to analyse multiple molecules, instead of repeatly calling function in `druglikeness.rulesfilter`, to save time

>>> rule = druglikeness.PC_rules(mols,n_jobs=4,detail=True)
>>> res = rule.CheckLipinskiRule()
>>> len(res)
760
>>> type(res)
<class 'list'>
>>> res[:3]
[{'MW': 256.26, 'logP': 2.83, 'nHD': 3, 'nHA': 3, 'Disposed': 'Accepted', 'nViolate': 0},
 {'MW': 288.25, 'logP': 2.79, 'nHD': 5, 'nHA': 5, 'Disposed': 'Accepted', 'nViolate': 0},
 {'MW': 182.17, 'logP': -3.59, 'nHD': 6, 'nHA': 6, 'Disposed': 'Accepted', 'nViolate': 1}]

Drug-likeness Filter
--------------------
The :mod:`structure_alert` module provides the tool to filter frequent hitters. This filter contains 11 endpoints

>>> from scopy.structure_alert import FilterWithSmarts

PAINS
~~~~~~
>>> res = FilterWithSmarts.Check_PAINS(mol)
>>> res
{'Disposed': 'Rejected', 'Endpoint': 'Pains'}

The function return a `dict`, the field :mod:`Disposed` represents compound state after filter applied (**Rejected** meant the compound rejected by filter, **Accepted** for accepted); :mod:`Endpoint` represents the which filter to be used.

Besides, the more specific information would be returned, if the :mod:`detail` has been set as :mod:`True` and the SMILES would be also returned if the :mod:`showSMILES` has been set as :mod:`True`.

>>> res = FilterWithSmarts.Check_PAINS(pains_mol, detail=True, showSMILES=True)
>>> res
{'SMILES': 'Cc1cc(O)cc(/N=C2/C=CC(=O)C(O)=C2)c1',
 'Disposed': 'Rejected',
 'MatchedAtoms': [((3, 2, 1, 0, 15, 16, 13, 12),)],
 'MatchedNames': ['Quinone_A'],
 'Endpoint': 'Pains'}

The result reveals the compound rejected by PAINS Filter, since the compound has the substructure named 'Quinone_A' which contained in PAINS Filter, more further, the No.3, No.2, No.1, No.0, No.15, No.16, No.13 and No.12 atom constructing this substructure.

Chelating
~~~~~~~~~
>>> res = FilterWithSmarts.Check_Chelating(pains_mol)
>>> res
{'Disposed': 'Accepted', 'Endpoint': 'Chelating'}

Promiscuity
~~~~~~~~~~~
>>> res = FilterWithSmarts.Check_BMS(pains_mol)
>>> res
{'Disposed': 'Accepted', 'Endpoint': 'BMS'}

Toxicity Filter
----------------
The :mod:`structure_alert.FilterWithSmarts` module also provides the tool to filter toxic compounds.

Acute Toxicity
~~~~~~~~~~~~~~
>>> res = FilterWithSmarts.Check_Genotoxic_Carcinogenicity_Mutagenicity(pains_mol)
>>> res
{'Disposed': 'Rejected', 'Endpoint': 'Genotoxic_Carcinogenicity_Mutagenicity'}

Samely, set :mod: `detail` to :mod: `True` to get specific infomation

>>> res = FilterWithSmarts.Check_Genotoxic_Carcinogenicity_Mutagenicity(pains_mol, detail=True)
>>> res
{'Disposed': 'Rejected',
 'MatchedAtoms': [((1, 0, 15, 16, 13), (12, 13, 15, 16, 0))],
 'MatchedNames': ['α, β-Unsaturated carbonyls'],
 'Endpoint': 'Genotoxic_Carcinogenicity_Mutagenicity'}

The molecule has matched the pattern twice

Environmental Toxicity
~~~~~~~~~~~~~~~~~~~~~~
Increasing attention to environmental impact of compounds in some regions.

>>> res = FilterWithSmarts.Check_NonBiodegradable(pains_mol, detail=True)
>>> res
{'Disposed': 'Rejected',
 'MatchedAtoms': [((0, 15, 16, 13),)],
 'MatchedNames': ['Ketone'],
 'Endpoint': 'NonBiodegradable'}

Multiprocessing
~~~~~~~~~~~~~~~
In reality, we trend to screen the compund library rather than sinlgle molecule. The :mod:`SmartsFilter` module provides the tool to screen multi-molecule under **Frequent Hitters Filter** and (or) **Toxicity Filter** 

>>> from scopy.structure_alert import SmartsFilter

>>> suppl = Chem.SDMolSupplier(os.path.join(ScoConfig.DemoDir,'760.sdf'))
>>> mols = (mol for mol in suppl if mol)

>>> F = SmartsFilter.Filter(mols, n_jobs=4, detail=True)
>>> res = F.Check_PAINS()
>>> type(res)
<class 'list'>

The function return a `list`
In above example, the PAINS Filter used to screen a library which contains 760 molecules under using four theardings

>>> res[0]
{'Disposed': 'Accepted', 'MatchedAtoms': ['-'], 'MatchedNames': ['-'], 'Endpoint': 'Pains'}
>>> res[207]
{'Disposed': 'Rejected', 'MatchedAtoms': [((7, 16, 15, 17, 18, 19, 20, 21, 14),)], 'MatchedNames': ['Mannich_A'], 'Endpoint': 'Pains'}

Chemical Space Analyser
-------------------------
To ensure obtaining a varity space of hitters, a Chemical space analysis of library is necessary before taking HTS. Chemical Space analysis could implement by calculating fingerprint (descriptor) and analysing framework (scaffold) of library.

Fingerprint Calculate
~~~~~~~~~~~~~~~~~~~~~
The :mod:`fingerprint` module provides the tool to compute fingerprints and (or) descriptor for chemical space analysis

EFG Fingerprints
"""""""""""""""""
Classification system termed “extended functional groups” (EFG), which are an extension of a set previously used by the CheckMol software, that covers in addition heterocyclic compound classes and periodic table groups. 

>>> from scopy.fingerprint import fingerprints
>>> fps = fingerprints.CalculateEFG(mols, useCount=False, n_jobs=4)
>>> fps.shape
(760, 583)
>>> fps.sum()
9473

In the above example, the calculated fingerprint is binary, beside that another type that using count to represent molecule(s)

>>> fps = fingerprints.CalculateEFG(mols, useCount=True, n_jobs=4)
>>> fps.shape
(760, 583)
>>> fps.sum()
58298

More details: `Salmina, Elena, Norbert Haider and Igor Tetko (2016)`_

.. _Salmina, Elena, Norbert Haider and Igor Tetko (2016):
	https://www.mdpi.com/1420-3049/21/1/1

IFG Fingerprint
""""""""""""""""
A new algorithm to identify all functional groups in organic molecules is presented.

>>> fps = fingerprints.CalculateIFG(mol, n_jobs=4)
>>> fps.shape
(760, 193)

Differ from other fingerprints, the dimension of IFG fingerprint may be variable with different library.

More details: `Peter Ertl (2017)`_.

.. _Peter Ertl (2017):
	https://jcheminf.springeropen.com/articles/10.1186/s13321-017-0225-z

Total 8 types of fingerprint are implemnted in :mod:`fingerprint`: MACCS, EFG, IFG, EState, Morgan, GhoseCrippen, Daylight and PubChem.

Framework Analyse
~~~~~~~~~~~~~~~~~
Beside counting frequency of each framework (scaffold), Scopy also supply a word cloud-like figure, called "Molecule Cloud".

The :mod:`visualize.mcloud` provide tool to implement framework analysis.

>>> from scopy.visualize import mcloud

Counting Frequency
""""""""""""""""""

>>> scount = mcloud.CountScaffold(mols)
>>> type(scount)
>>> dict
>>> scount['c1ccccc1']
>>> 68

The function return a `dict` whose `keys` is the scaffold in SMILES format and `values` is the corresponding frequnecy

Molecule Cloud
"""""""""""""""
.. note::
	This module should run under a Java environment and the script retrived from `Peter Ertl`_

>>> scaffolds = os.path.join(ScoConfig.DemoDir, 'scaffolds.txt')
>>> mcloud.ShowMcloud(file=scaffolds, number=200, savedir='./mcloud.png')

.. figure:: /image/mcloud.png
	:width: 400px
	:align: center
	
	Molecule cloud, more frequent molecule (scaffold) appear, the more bigger and more forward layer got.


Screening Visualizer
--------------------
The :mod:`visualize` module provides the tool to visualize PC properties, PC-drived rules, substructures, fingerprints and molecular scaffolds (see `Molecule Cloud`_).

PC Visualizer
~~~~~~~~~~~~~
The :mod:`visualize.pc_depict` module can depict basic properties distribution of molecule(s) and position molecular values within the selected filter range.

Proprty Matrix
""""""""""""""
The proprty matrix can intuitively show the compounds' distribution in Two-Dimension space, and diagonal of the matrix is the displot of property

>>> from scopy.visualize import pc_depict
>>> fig = pc_depict.prop_matrix()
>>> fig
<Figure size 1567x989 with 36 Axes>

.. figure:: /image/760_matrix.png
	:width: 400px
	:align: center

	The matrix of logP, TPSA, MW, nRot, nHD and nHA

Default properties of matrix are logP, TPSA, MW, nRot, nHD and nHA. The user could customize proerties to be shown through parament `items`

>>> fig = pc_depict.prop_matrix(mols, n_jobs=4, items=['MW', 'Vol', 'Dense'])

.. figure:: /image/760_matrix_2.png
	:width: 400px
	:align: center

	The matrix of MW, Vol and Dense

Basic Property Radar
""""""""""""""""""""
 A radar plot positionning compound's values within the selected filter ranges (pale blue and red). By default, the `drug-like soft`_ filter ranges are visualized.

.. note::
	The property "Number of Charged Groups" in `drug-like soft`_ has not been implemented

>>> fig = pc_depict.rule_radar(mol)
>>> fig
<Figure size 640x480 with 1 Axes>

.. figure:: /image/mol_basci_rule.png
	:width: 400px
	:align: center
	
	A radar plot of drug-like soft

.. _`drug-like soft`: http://fafdrugs4.mti.univ-paris-diderot.fr/filters.html

Fragment visualizer
~~~~~~~~~~~~~~~~~~~
The :mod:`visualize.highlight` module can flag the subtructure related to some specific endpoint

>>> from scopy.visualize import highlight

>>> fig = highlight.HighlightAtoms(mol,highlightAtoms=[3, 2, 1, 0, 15, 16, 13, 12]) 
>>> type(fig)
IPython.core.display.SVG

The atoms highlighted retrieved from module :mod:`structure_alert.FilterWithSmarts` and the function return a SVG object, you may should save it manually.

>>> with open(r'PAINS.svg', 'w') as f_obj:
		f_obj.write(fig.data)
	f_obj.close()

.. figure:: /image/PAINS.svg
	:width: 400px
	:align: center

	A molecule with highlighted substructure in red


 









