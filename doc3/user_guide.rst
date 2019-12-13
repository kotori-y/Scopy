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

Install from source
~~~~~~~~~~~~~~~~~~~
>>> git clone git@github.com:kotori-y/Scopy.git && cd scopy
>>> [sudo] python setup.py install

PyPI
~~~~
>>> 

Conda
~~~~~
>>> 

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

You could calculate different properties respectively

>>> MW = molproperty.CalculateMolWeight(mol)
>>> MW
229.24
>>> logP = molproperty.CalculateLogP(mol)
>>> logP
2.35
>>> nHD = molproperty.CalculateNumHDonors(mol)
>>> nHD
2

You could also calculate total(43) properties at once time

>>> props = molproperty.GetProperties(mol)
>>> props
{'MW': 229.24, 'Vol': 235.2, 'Dense': 0.97, 'fChar': 0, 'nBond': 18, 'nAtom': 28, 'nHD': 2, 'nHA': 4, 'nHB': 6, 'nHet': 4, 'nStero': 0, 'nHev': 17, 'nRot': 1, 'nRig': 14, 'nRing': 2, 'logP': 2.35, 'logD': 0.8670776309515202, 'pKa': -5.931602224875785, 'logSw': -2.95, 'ab': 'base', 'MR': 64.8, 'tPSA': 69.89, 'AP': 0.35, 'HetRatio': 0.31, 'Fsp3': 0.08, 'MaxRing': 6, 'QEDmean': 0.73, 'QEDmax': 0.7, 'QEDnone': 0.79, 'SAscore': 2.96, 'NPscore': 0.49, 'nSingle': 8, 'nDouble': 4, 'nTriple': 0, 'nC': 13, 'nB': 0, 'nF': 0, 'nCl': 0, 'nBr': 0, 'nI': 0, 'nP': 0, 'nS': 0, 'nO': 3, 'nN': 1}

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

In above example, the compound do not violate any property limited in Lipinski Rule and its status is 'Accepted'.

Besides, the specific value of each propety would be returned if the :mod:`detail` has been set as :mod:`True`.

>>> res = rulesfilter.CheckLipinskiRule(mol, detail=True)
>>> res
{'MW': 229.24, 'logP': 2.35, 'nHD': 2, 'nHA': 4, 'Disposed': 'Accepted', 'nViolate': 0}

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

Structure Alert Filter
----------------------
The :mod:`structure_alert` module provides the tool to search for the presence of toxicophores and flag unwanted reactive chemical groups, where both toxicophores and unwanted reactive chemical groups are encoded by SMARTS.

Screening a single molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`FilterWithSmarts` module provides the tool to screen a molecule. 

>>> from scopy.structure_alert import FilterWithSmarts
>>> mol = Chem.MolFromMolFile(os.path.join(ScoConfig.DemoDir,'PAINS.sdf'))
>>> mol
<rdkit.Chem.rdchem.Mol object at 0x000001E75CE60580>

In here, the PAINS Filter be used to screen a molecule

>>> res = FilterWithSmarts.Check_PAINS(mol)
>>> res
{'Disposed': 'Rejected', 'Endpoint': 'Pains'}

The function return a `dict`, the field :mod:`Disposed` represents compound state after filter applied (**Rejected** meant the compound rejected by filter, **Accepted** for accepted); :mod:`Endpoint` represents the which filter to be used.

Besides, the more specific information would be returned, if the :mod:`detail` has been set as :mod:`True` (defaults to :mod:`False`)

>>> res = FilterWithSmarts.Check_PAINS(mol, detail=True)
>>> res
{'Disposed': 'Rejected', 'MatchedAtoms': [((3, 2, 1, 0, 15, 16, 13, 12),)], 
 'MatchedNames': ['Quinone_A'], 'Endpoint': 'Pains'}

The result reveals the compound rejected by PAINS Filter, since the compound has the substructure named 'Quinone_A' which contained in PAINS Filter, more further, the No.3, No.2, No.1, No.0, No.15, No.16, No.13 and No.12 atom constructing this substructure.

Screening multi-molecule
~~~~~~~~~~~~~~~~~~~~~~~~
In reality, we trend to screen the compund library rather than sinlgle molecule. The :mod:`SmartsFilter` module provides the tool to screen multi-molecule 

>>> from scopy.structure_alert import SmartsFilter

>>> suppl = Chem.SDMolSupplier(os.path.join(ScoConfig.DemoDir,'760.sdf'))
>>> mols = [mol for mol in suppl if mol]

>>> F = SmartsFilter.Filter(mols, n_jobs=4, detail=True)
>>> res = F.Check_PAINS()

In above example, the PAINS Filter used to screen a library which contains 760 molecules under using four theardings

The function return a `list`

>>> res[0]
{'Disposed': 'Accepted', 'MatchedAtoms': ['-'], 'MatchedNames': ['-'], 'Endpoint': 'Pains'}
>>> res[207]
{'Disposed': 'Rejected', 'MatchedAtoms': [((7, 16, 15, 17, 18, 19, 20, 21, 14),)], 'MatchedNames': ['Mannich_A'], 'Endpoint': 'Pains'}

Fingerprint Calculator
----------------------
The :mod:`fingerprint` module provides the tool to compute fingerprints retrieved from fragments.

EFG Fingerprint
~~~~~~~~~~~~~~~
Classification system termed “extended functional groups” (EFG), which are an extension of a set previously used by the CheckMol software, that covers in addition heterocyclic compound classes and periodic table groups. 

>>> from scopy.fingerprint import fingerprints
>>> fps = fingerprints.CalculateEFG(mol, useCount=False, n_jobs=4)
>>> fps.shape
(760, 583)
>>> fps.sum()
9473

In the above example, the calculated fingerprint is binary, beside that another type that using count to represent molecule(s)

>>> fps = fingerprints.CalculateEFG(mol, useCount=True, n_jobs=4)
>>> fps.shape
(760, 583)
>>> fps.sum()
58298

More details: `Salmina, Elena, Norbert Haider and Igor Tetko (2016)`_

.. _Salmina, Elena, Norbert Haider and Igor Tetko (2016):
	https://www.mdpi.com/1420-3049/21/1/1

IFG Fingerprint
~~~~~~~~~~~~~~~
A new algorithm to identify all functional groups in organic molecules is presented.

>>> fps = fingerprints.CalculateIFG(mol, n_jobs=4)
>>> fps.shape
(760, 193)

Differ from other fingerprints, the dimension of IFG fingerprint may be variable with different library.

More details: `Peter Ertl (2017)`_.

.. _Peter Ertl (2017):
	https://jcheminf.springeropen.com/articles/10.1186/s13321-017-0225-z


Total 8 types of fingerprint are implemnted in :mod:`fingerprint`: MACCS, EFG, IFG, EState, Morgan, GhoseCrippen, Daylight and PubChem.

Screening Visualizer
--------------------
The :mod:`visualize` module provides the tool to visualize PC properties, PC-drived rules, substructures, fingerprints and molecular scaffolds.

PC Visualizer
~~~~~~~~~~~~~
The :mod:`visualize.pc_depict` module can depict basic properties distribution of molecule(s) and position molecular values within the selected filter range.

Proprty Matrix
""""""""""""""
The proprty matrix can intuitively show the compounds' distribution in Two-Dimension space, and diagonal of the matrix is the displot of property

>>> from scopy.visualize import pc_decipt
>>> pc = pc_depict.PropVisual()
>>> fig = pc.prop_matrix(mols, n_jobs=4)
>>> fig
<Figure size 1567x989 with 36 Axes>

.. figure:: /image/760_matrix.png
	:width: 600px
	:align: left
	The matrix of logP, TPSA, MW, nRot, nHD and nHA

Default properties of matrix are logP, TPSA, MW, nRot, nHD and nHA. The user could customize proerties to be shown through parament `items`

>>> fig = pc.prop_matrix(mols, n_jobs=4, items=['MW', 'Vol', 'Dense'])

.. figure:: /image/760_matrix_2.png
	:width: 600px
	:align: left
	The matrix of MW, Vol and Dense

Basic Property Radar
""""""""""""""""""""
 A radar plot positionning compound's values within the selected filter ranges (pale blue and red). By default, the `drug-like soft`_ filter ranges are visualized.

.. note::
	The property "Number of Charged Groups" in `drug-like soft`_ has not been implemented

>>> fig = pc.rule_radar(mol)
>>> fig
<Figure size 640x480 with 1 Axes>

.. figure:: /image/mol_basci_rule.png
	:width: 600px
	:align: left
	A radar plot of drug-like soft

.. _`drug-like soft`: http://fafdrugs4.mti.univ-paris-diderot.fr/filters.html

	
Fragment visualizer
~~~~~~~~~~~~~~~~~~~


Fingerprint visualizer
~~~~~~~~~~~~~~~~~~~~~

Molcule Cloud
~~~~~~~~~~~~~

Molecular Pretreater
-------------------





.. figure:: /image/CPI.png
	:width: 400px
	:align: center
	
	The calculation process for chemical-protein interaction descriptors.