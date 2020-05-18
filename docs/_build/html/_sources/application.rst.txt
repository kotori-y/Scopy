Application
===========
In order to obtain a better understanding of Scopy and its application, we introduce an examples to show a regular screening process. The data in this application can be obtain at https://github.com/kotori-y/Scopy_Data.

.. figure:: /image/frame.png
	:width: 600px
	:align: center 
	
	The overview of Scopy python package.

Application Screening of Small Molecule Libarary
--------------------------------------------------
We had collected 408,598 molecules from Zelin database. At first, we took a framework analysis of library through molecular cloud, and then 3 drug-likeness filters (Lipinski, Pfizer and GSK), 3 toxicity endpoint filters (Potential_Electrophilic, Skin_Sensitization and LD50_oral Filter) and 6 frequent hit compound filters (AlphaScreen_FHs, AlphaScreen_GST_FHs, AlphaScreen_HIS_FHs, Chelating, BMS and PAINS Filter) were used for screening library.

.. code-block:: python
	:linenos:

	from rdkit import Chem
	import pandas as pd #This package should be installed
	import matplotlib.pyplot as plt
	from scopy.Scopretreat import pretreat
	from scopy.ScoVisualize import mcloud, highlight, pc_depict
	from scopy.ScoFH import FHfilter
	from scopy.ScoTox import Toxfilter
	from scopy.ScoDruglikeness import PC_properties, PC_rules
	   
	#=========================================================================
	# Read library
	#=========================================================================
	data = pd.read_csv('Zelin.csv')
	mols = [Chem.MolFromSmiles(smi) for smi in data.mol.values]
	
	#=========================================================================
	# Standardize molecules
	#=========================================================================
	mols = [pretreat.StandardMol(mol) for mol in mols]

	#=========================================================================
	# The instantiate step and the use of 20 processers.
	#=========================================================================
	props = PC_properties(mols, n_jobs=20)
	rules = PC_rules(mols, n_jobs=20, detail=True)
	fh = FHfilter(mols, n_jobs=20, detail=True)
	tox = Toxfilter(mols, n_jobs=20, detail=True)

	#=========================================================================
	# Framework analyse
	#=========================================================================
	scount = mcloud.CountScaffold(mols,stype='Murcko')
	scount = pd.DataFrame(scount,index=['Frequency']).T
	#scount.to_csv('scount.txt',sep='\t',header=None)

	#=========================================================================
	# PC Properties.
	#=========================================================================
	prop_res = props.GetProperties(items=['MW','logP','nHA','nHD','TPSA']) #5 propperties are choossen
	prop_res = pd.DataFrame(prop_res)

	#=========================================================================
	# Drug-likeness rule
	#=========================================================================
	gsk = pd.DataFrame(rules.CheckGSKRule())
	pizer = pd.DataFrame(rules.CheckPfizerRule())
	lipiski = pd.DataFrame(rules.CheckLipinskiRule())
	   
	summary_druglike = pd.DataFrame({'SMILES':data.mol.values,
	                                 'Pizer':pizer.Disposed,
	                                 'GSK':gsk.Disposed})

	summary_druglike['Rejected_Num'] = (summary_druglike=='Rejected').sum(axis=1)
	summary_2_druglike = pd.DataFrame((summary_druglike.iloc[:,1:]=='Rejected').sum(axis=0))
	summary_2_druglike.columns = ['Rejected']
	summary_2_druglike['Accepted'] = len(summary_druglike)-summary_2_druglike.Rejected.values

	#=========================================================================
	# Toxicity
	#=========================================================================
	ele = pd.DataFrame(tox.Check_Potential_Electrophilic())
	skin = pd.DataFrame(tox.Check_Skin_Sensitization())
	ld_50 = pd.DataFrame(tox.Check_LD50_Oral())
	   
	summary_tox = pd.DataFrame({'SMILES':data.mol.values,
	                            'Potential_Electrophilic':ele.Disposed,
	                            'Skin_Sensitization':skin.Disposed,
	                            'LD50_Oral':ld_50.Disposed,
	                            })
	   
	summary_tox['Rejected_Num'] = (summary_tox=='Rejected').sum(axis=1)
	summary_2_tox = pd.DataFrame((summary_tox.iloc[:,1:]=='Rejected').sum(axis=0))
	summary_2_tox.columns = ['Rejected']
	summary_2_tox['Accepted'] = len(summary_tox)-summary_2_tox.Rejected.values
	   
	#=========================================================================
	# Frequent Hitters
	#=========================================================================
	alapha_fh = pd.DataFrame(fh.Check_AlphaScreen_FHs())
	gst = pd.DataFrame(fh.Check_AlphaScreen_GST_FHs())
	his = pd.DataFrame(fh.Check_AlphaScreen_HIS_FHs())
	che = pd.DataFrame(fh.Check_Chelating())
	bms = pd.DataFrame(fh.Check_BMS())
	pains = pd.DataFrame(fh.Check_PAINS())


	summary_fh = pd.DataFrame({'SMILES':data.mol.values,
	                           'AlphaScreen_FHs':alapha_fh.Disposed,
	                           'AlphaScreen_GST_FHs':gst.Disposed,
	                           'AlphaScreen_HIS_FHs':his.Disposed,
	                           'Chelating':che.Disposed,
	                           'BMS':bms.Disposed,
	                           'PAINS':bms.Disposed})

	summary_fh['Rejected_Num'] = (summary_fh=='Rejected').sum(axis=1)
	summary_2_fh = pd.DataFrame((summary_fh.iloc[:,1:]=='Rejected').sum(axis=0))
	summary_2_fh.columns = ['Rejected']
	summary_2_fh['Accepted'] = len(summary_fh)-summary_2_fh.Rejected.values

The summary of the Murcko scaffolds of the database.

>>> scount
                                             Frequency
O=C(Nc1cn[nH]c1)c1cn[nH]c1                          26
O=C(Nc1ccccc1)c1nc2ncccn2n1                        151
c1csc(-c2ccc3ccccc3n2)c1                           117
c1ccccc1                                         13723
O=C(c1cc[nH]n1)N1CCN(c2ncccn2)CC1                    8
...                                                ...
[NH2+]=C(Nc1ncccn1)N(c1ccccc1)S(=O)(=O)c1ccccc1      1
O=C1C(C=[NH+]Cc2ccccc2)=Cc2ccccc21                   1
N=C(NC(=O)c1ccccc1)Nc1nc2ccccc2s1                    1
O=S(=O)(c1cccc2nsnc12)n1c[nH+]c2ccccc21              1
O=C(C=Cc1cccs1)Nc1nnc[nH]1                           1

The cloud gram of the top 200 Murcko scaffolds appeared in the database.

>>> mcloud.ShowMcloud(r"scount.txt",skip=1,number=200) #The skip parameter is used to skip the most frequent skeleton (benzene ring)

.. figure:: /image/application/mcloud.png
	:width: 600px
	:align: center

The values of five physicochemical properties of the database, including MW, logP, nHA, nHD and TPSA.

>>> prop_res
            MW  logP  nHA  nHD    TPSA
0       290.15  0.78    6    2  107.83
1       355.09  0.95    9    1  124.78
2       359.97  3.65    4    0   53.02
3       281.00  3.54    2    1   46.17
4       272.14  0.17    6    0   67.15
...        ...   ...  ...  ...     ...
408592  433.05  4.95    3    3   59.37
408593  273.06  3.20    5    1   55.98
408594  346.17  1.30    2    2   36.78
408595  367.21  4.62    4    0   55.73
408596  433.05  4.95    3    3   59.37

The feature-feature related scatter diagram of five features of the database, including MW, logP, nHA, nHD and TPSA.

>>> fig_matrix = pc_depict.prop_matrix(mols, items=['MW','logP','nHA','nHD','TPSA'], n_jobs=20)

.. figure:: /image/application/matrix.png
	:width: 500px
	:align: center

The radar plot of molecule 0 (CCn1cc(C(=O)Nc2cn(CC)nc2C(N)=O)c(C)n1).

>>> prop_kws = {'MW':[0,500],'logP':[None,5],'nHD':[0,5],
	            'nHA':[0,10],'TPSA':[0,140]} #Linpiski's Rule and TPSA from Veber's Rule
>>> fig_radar = pc_depict.rule_radar(mol, prop_kws)

.. figure:: /image/application/radar.png
	:width: 400px
	:align: center

The summary of screening results.

>>> summary_2_druglike
     Filter  Rejected  Accepted
0  Lipinski      8448    400149
1     Pizer     73193    335404
2       GSK    200192    208405
>>>#
>>> summary_2_tox
                    Filter  Rejected  Accepted
0  Potential_Electrophilic    111946    296651
1       Skin_Sensitization    310897     97700
2                LD50_Oral      9466    399131
>>>#
>>>summary_2_fh
                Filter  Rejected  Accepted
0      AlphaScreen_FHs         3    408594
1  AlphaScreen_GST_FHs      3500    405097
2  AlphaScreen_HIS_FHs      2954    405643
3            Chelating      6757    401840
4                  BMS     92553    316044
5                PAINS     19375    389222

And the visualization of screening result.

.. code-block:: python
	:linenos:

	f,axes = plt.subplots(1,3,figsize=(5*3,9))

	labels = ['six-filter','five-filter','four-filter','three-filter','two-filter','one-filter','zero-filter']

	res_druglike = summary_druglike.Rejected_Num.value_counts().sort_index(ascending=False)
	res_tox = summary_tox.Rejected_Num.value_counts().sort_index(ascending=False)
	res_fh = summary_fh.Rejected_Num.value_counts().sort_index(ascending=False)

	axes[0].pie(res_druglike.values,explode=(0,0,0,0.1),
	            labels=labels[-4:],autopct='%1.1f%%',
	            shadow=False,startangle=150)
	axes[0].set_title('Drug-likeness Rule Filter')

	axes[1].pie(res_tox.values,explode=(0,0,0,0.1),
	            labels=labels[-4:],autopct='%1.1f%%',
	            shadow=False,startangle=150)
	axes[1].set_title('Toxicity Filter')

	axes[2].pie(res_fh.values,explode=(0,0,0,0,0.1),
	            labels=labels[-5:],autopct='%1.1f%%',
	            shadow=False,startangle=150)
	axes[2].set_title('FH Filter')

	plt.show()   
    

.. figure:: /image/application/small_mol_ratio.png
	:width: 600px
	:align: center

From the figure, it is shown that 37.6% of small molecules failed the drug-likeness screening, most of which were filtered by the GSK rule. One of the main reasons is that the MW limitation of the GSK rule is 400, which is strict for commercial molecule database. About 30% of the molecules failed the toxicity filters, such as "skin_sensitvity endpoint", which can be avoided by changing the route of administration; in the end, nearly three-quarters of the molecules passed the frequent hit compound filter.

Take molecule_21 (O=C(N/N=C/c1cc(O)c(O)cc1)[C@@H]1[C@H](c2ccc(C(C)(C)C)cc2)C1) as an example to show the functional group marker gram.

>>> from scopy.structure_alert import FilterWithSmarts
>>>	
>>> mol = Chem.MolFromSmiles('O=C(N/N=C/c1cc(O)c(O)cc1)[C@@H]1[C@H](c2ccc(C(C)(C)C)cc2)C1')
>>> mol = pretreat.StandardMol(mol)
>>> FilterWithSmarts.Check_PAINS(mol,detail=True)
{'Disposed': 'Rejected',
 'MatchedAtoms': [((9, 7, 6, 5, 12, 11, 4, 3, 2, 10),),
  ((5, 12, 11, 9, 7, 6, 8, 10),)],
 'MatchedNames': ['Hzone_phenol_B', 'Catechol_A'],
 'Endpoint': 'Pains'}

>>> highlight.HighlightAtoms(mol,highlightAtoms=[9, 7, 6, 5, 12, 11, 4, 3, 2, 10])

.. figure:: /image/application/highlight_1.svg
	:width: 500px
	:align: center

>>> highlight.HighlightAtoms(mol,highlightAtoms=[5, 12, 11, 9, 7, 6, 8, 10])

.. figure:: /image/application/highlight_2.svg
	:width: 500px
	:align: center