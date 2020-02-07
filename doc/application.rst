Application
===========
The Scopy Python package can filter compounds library in three aspects: drug-likeness, frequent hit and potential toxicity. Poor drug-likeness make it difficult for compounds to become drugs, and frequent hitters, including promiscuous ligand and flase positive compounds, would stain our hitters. Besides, the toxic compound should be removed from our lirary before taking HTS. We will introduce **X** examples of Scopy's applications including ....... to show some plan of pre-screening in different situation before formal HTS screening. You could download these scripts via: xxxx, and data: xxxxx

.. figure:: /image/frame.png
	:width: 400px
	:align: center 
	
	The overview of Scopy python package. Scopy can filter compounds with three filters: drug-likeness, frequent hitters and toxicity filters. Besides, the Chemical Space Analyser provide users with a general understanding of the chemical space of screening libraries. The Visualizer allows users to visually analyze the results.


Application 1. Screening Withdrawn Drugs  
-----------------------------------------
Drug withdrawn due to adverse reactions caused by drug toxicity is a huge loss in the drug development process, often accompanied by long legal proceedings and huge compensation. If we colud filter these potential toxic compounds, the withdrawal would be avoided. 
We had obtained **222** withdrawn drugs from `DrugBank`_ and used filter to screen them, and the toxicity filter has been choosen. We have choosen broad toxicity endpoint (Potential_Electrophilic Filter), acute toxicity endpoint (including LD50_oral, Genotoxic_Carcinogenicity_Mutagenicity and NonGenotoxic_Carcinogenicity Filter) and comprehensive toxicity endpoint (NTD, SureChEMBL and Toxicophores filter)

.. _`DrugBank`: https://www.drugbank.ca/

.. code-block:: python
	:linenos:

	import pandas as pd #This package should be installed
	from rdkit import Chem
	from scopy.structure_alert.SmartsFilter import Filter


	data = pd.read_csv(r"withdrawn_drug.csv")    
	mols = [Chem.MolFromSmiles(smi) for smi in data.mol.values]

	screener= Filter(mols,n_jobs=4,detail=True)#Instantiate

	#=========================================================================
	# Broad Toxicity
	#=========================================================================
	ele_res = pd.DataFrame(screener.Check_Potential_Electrophilic()) #Potential_Electrophilic

	#=========================================================================
	# Acute Toxicity
	#=========================================================================
	ld50_res = pd.DataFrame(screener.Check_LD50_Oral()) #LD50_Oral
	gene_res = pd.DataFrame(screener.Check_Genotoxic_Carcinogenicity_Mutagenicity()) #Genotoxic_Carcinogenicity_Mutagenicity
	nogene_res = pd.DataFrame(screener.Check_NonGenotoxic_Carcinogenicity()) #NonGenotoxic_Carcinogenicity

	#=========================================================================
	# Comprhenesive
	#=========================================================================
	ntd_res = pd.DataFrame(screener.Check_NTD()) #NTD
	toxci_res = pd.DataFrame(screener.Check_NonGenotoxic_Carcinogenicity()) #Toxicophores
	chemble_res = pd.DataFrame(screener.Check_SureChEMBL()) #NonGenotoxic_Carcinogenicity

	#=========================================================================
	# Get Summary Result
	#=========================================================================
	summary = pd.DataFrame({'SMILES': data.mol.values,
	                        'Potential_Electrophilic': ele_res.Disposed,
	                        'LD50_Oral': ld50_res.Disposed,
	                        'Genotoxic_Carcinogenicity_Mutagenicity': gene_res.Disposed,
	                        'NonGenotoxic_Carcinogenicity': nogene_res.Disposed,
	                        'NTD': ntd_res.Disposed,
	                        'Toxicophores': toxci_res.Disposed,
	                        'SureChEMBL': chemble_res.Disposed})
	summary['Rejected_Num'] = (summary == 'Rejected').sum(axis=1)

	summary_2 = pd.DataFrame((summary.iloc[:,1:-1] == 'Rejected').sum(axis=0), columns=['Rejected'])
	summary_2['Accepted'] = 222 - summary_2.Rejected.values

>>> summary
                                                SMILES Genotoxic_Carcinogenicity_Mutagenicity  ... Toxicophores Rejected_Num
0    OC[C@@H](C(=O)N([C@H](C(=O)N[C@@H](C(=O)N[C@H]...                               Rejected  ...     Rejected            5
1                      OC(=O)[C@H](Cc1c[nH]c2c1cccc2)N                               Accepted  ...     Rejected            1
2                                   CCC(C#C)(/C=C\Cl)O                               Rejected  ...     Rejected            5
3    O=C1NC(=O)C(S1)Cc1ccc(cc1)OCC1(C)CCc2c(O1)c(C)...                               Accepted  ...     Rejected            3
4        COc1ccc(cc1)C(=C(c1ccc(cc1)OC)Cl)c1ccc(cc1)OC                               Accepted  ...     Rejected            4
..                                                 ...                                    ...  ...          ...          ...
217                            CCCCOc1ccc(cc1)CC(=O)NO                               Rejected  ...     Rejected            4
218     OC[C@H]1O[C@H](C[C@@H]1O)n1cc(CC)c(=O)[nH]c1=O                               Accepted  ...     Accepted            0
219                   OC(=O)COc1nn(c2c1cccc2)Cc1ccccc1                               Accepted  ...     Rejected            1
220                   CCN(CCOc1ccc2c(c1)sc(n2)N(C)C)CC                               Rejected  ...     Rejected            2
221      Oc1cc2O[C@H](c3ccc(c(c3)O)O)[C@H](Cc2c(c1)O)O                               Accepted  ...     Rejected            3

>>> summary_2
                                   Filter  Rejected  Accepted
0  Genotoxic_Carcinogenicity_Mutagenicity        94       128
1                               LD50_Oral        13       209
2            NonGenotoxic_Carcinogenicity        67       155
3                                     NTD        99       123
4                 Potential_Electrophilic        46       176
5                              SureChEMBL        39       183
6                            Toxicophores       157        65

.. code-block:: python
	:linenos:

	import matplotlib.pyplot as plt

	f,ax = plt.subplots()

	labels = ['six-filter','five-filter','four-filter','three-filter','two-filter','one-filter','zero-filter']
	res = summary.Rejected_Num.value_counts().sort_index(ascending=False)
	sizes = res.values
	explode = (0,0,0,0,0,0,0.1)

	ax.pie(sizes,explode=explode,labels=labels,autopct='%1.1f%%',shadow=False,startangle=150)
	plt.show()  

.. figure:: /image/withdrawn_ratio.png
	:width: 400px
	:align: center 
	
	The ratio of the withdrawn drugs detected by specific number filter.

From the above pie chart, we can see 89.2% of withdrawn drugs has been detected at least one filter, and the `summary_2` presents most drugs could be filtered by Toxicophores endpoint.

Next, we take Tolrestat (COc1ccc2c(c1C(F)(F)F)cccc2C(=S)N(CC(=O)O)C) and Toxicophores as an example to show which substructure has been regarded as 'unwanted' groups

>>> from rdkit import Chem
>>> from scopy.structure_alert import FilterWithSmarts

>>> mol = Chem.MolFromSmiles('COc1ccc2c(c1C(F)(F)F)cccc2C(=S)N(CC(=O)O)C')
>>> FilterWithSmarts.Check_Toxicophores(mol, detail=True)
{'Disposed': 'Rejected',
 'MatchedAtoms': [((15, 16, 17, 18, 19, 23),), ((16, 17),), ((17, 16),)],
 'MatchedNames': ['thioamide', 'thiocarbonyl_aromatic', 'thioketone'],
 'Endpoint': 'Toxicophores'}

The result shows that the Tolrestat contain three substructures (thioamide, thiocarbonyl_aromatic, thioketone) which be regarded as 'Toxicophores' groups, and the 'thioamide' substructure was consisted of No.15, 16, 17, 18, 19 and 23 atoms

>>> from scopy.visualize import highlight
>>> highlight.HighlightAtoms(mol, highlightAtoms=(15, 16, 17, 18, 19, 23))

.. figure:: /image/ex_withdrawn.svg
	:width: 400px
	:align: center



Application 2. Screening Natural Products  
-------------------------------------------
Natural products are rich in structure, so if we want to screen natural product libraries, a chemical spatial analysis is necessary. It is then screened for drug-like properties and toxicity. 
We collected 20,015 natural products to form a screening library. First, the molecular library was used for scaffold analysis of the screening library, then four drug-likeness rules (Lipinski Rule, bRo5 Rule, Xu's Rule and Oral Macrocycle Rule) were used for drug-like analysis, and tree comprehensive toxicity filters (NTD, SureChEMBL and Toxicophores filter) were used for toxicity analysis. Finally calculate six commonly used physical and chemical properties (MW, logP, nHA, nHD, TPSA and nRot) for subsequent analysis.

.. code-block:: python
	:linenos:

	from rdkit import Chem
	import pandas as pd #This package should be installed
	from scopy.structure_alert.SmartsFilter import Filter
	from scopy.druglikeness.druglikeness import PC_properties, PC_rules
	from scopy.visualize import mcloud

	data = pd.read_csv(r'Natural_product-NPS20015.csv')
	mols = [Chem.MolFromSmiles(smi) for smi in data.Canonical_SMILES.values]

	#=========================================================================
	# Statistics of Murcko's framework frequency
	#=========================================================================
	scount = mcloud.CountScaffold(mols,stype='Murcko')
	scount = pd.DataFrame(scount,index=['Frequency']).T
	#scount.to_csv('scount.txt',sep='\t',header=None)
		
	props = PC_properties(mols,n_jobs=-1) #Calculate physical and chemical properties
	rules = PC_rules(mols,n_jobs=-1,detail=True) #Screen drug-likeness rules
	screener = Filter(mols,n_jobs=4,detail=True) #Screen Toxicity

	#=========================================================================
	# Drug-likeness rule
	#=========================================================================
	bRo5_res = pd.DataFrame(rules.CheckBeyondRo5())
	Macro_res = pd.DataFrame(rules.CheckOralMacrocycles())
	lipiski_res = pd.DataFrame(rules.CheckLipinskiRule())
	xu_res = pd.DataFrame(rules.CheckLipinskiRule())

	#=========================================================================
	# Comprhenesive Toxicity
	#=========================================================================
	ntd_res = pd.DataFrame(screener.Check_NTD())
	toxci_res = pd.DataFrame(screener.Check_Toxicophores())
	chemble_res = pd.DataFrame(screener.Check_SureChEMBL())

	#=========================================================================
	# PC properties
	#=========================================================================
	MW = props.CalculateMolWeight()
	logP = props.CalculateLogP()
	nHA = props.CalculateNumHAcceptors()
	nHD = props.CalculateNumHDonors()
	tPSA = props.CalculateTPSA()
	nRot = props.CalculateNumRotatableBonds()

	#=========================================================================
	# Summary
	#=========================================================================

	summary = pd.DataFrame({'mol': data.Canonical_SMILES.values,
	                        'bRo5': bRo5_res.Disposed,
		                    'Macro': Macro_res.Disposed,
	                        'Lipinski': lipiski_res.Disposed,
	                        'Xu': Xu_res.Disposed,
	                        'NTD': ntd_res.Disposed,
	                        'Toxicophores': toxci_res.Disposed,
	                        'SureChemble': chemble_res.Disposed,
	                        'MW': MW,
	                        'logP':log
	                        'nHA': nHA,
	                        'nHD':nHD,
	                        'TPSA':TPSA,
	                        'nRot':nRot})
	summary['Rejected_Num'] = (summary == 'Rejected').sum(axis=1)
	summary_2 = pd.DataFrame((summary.iloc[:,1:-7] == 'Rejected').sum(axis=0), columns=['Rejected'])
	summary_2['Accepted'] = 20015 - summary_2.Rejected.values

>>> scount
                                                    Frequency
c1ccccc1                                                 1227
C1=C2C3CCCCC3CCC2C2CCC3CCCCC3C2C1                         376
O=c1cc(-c2ccccc2)oc2ccccc12                               333
C1=CCCCC1                                                 209
C1CCC2C(C1)CCC1C3CCCC3CCC21                               203
O=C1CC(c2ccccc2)Oc2ccccc21                                175
...                                                       ...                          
C=C1COC2COC(OCC34CC5CCCC5C5CC3CCC54)CC12                    1
O=C1CC2CCCC3N4CCCC1C23CCC4                                  1
O=C1CC2C(CCC3C4C=C(CCCCOC5CCCCO5)OC4CC32)C2CCC(...          1
O=C1CC2OCC3CCC4C(CCC5C(c6ccoc6)OC(=O)C6OC654)C3...          1
O=C1OCC2OC(OC(=O)c3ccccc3-c3ccccc31)C1CC2Oc2ccc...          1
C=C1CCC2C=COCC12                                            1

Next we use molecular cloud to visualize the results intuitively

>>> mcloud.ShowMcloud(r"scount.txt",skip=1,number=200) #The skip parameter is used to skip the most frequent skeleton (benzene ring)

.. figure:: /image/natural_mcloud.png
	:width: 400px
	:align: center

>>> summary
                                                     mol      bRo5     Macro  Lipinski  ... nHD    tPSA nRot Rejected_Num
0                Oc1cc2[OH+][C-](c3ccccc3)C(=Cc2c(c1)O)O  Accepted  Accepted  Accepted  ...   3   73.49    1            2
1       OC1=C[C-]2[OH+]C(=C(C=C2C(=C1)O)O)c1ccc(c(c1)O)O  Accepted  Rejected  Accepted  ...   5  113.95    1            4
2                                   OCC(C(C(C(CO)O)O)O)O  Rejected  Rejected  Accepted  ...   6  121.38    5            2
3      Oc1cc(O)c2c(c1)OC(C(C2c1c(O)cc(c2c1OC(C(C2)O)c...  Rejected  Rejected  Rejected  ...  10  220.76    3            6
4      COc1c(O)cc(cc1O)C1Oc2c(CC1O)c(O)cc(c2C1C(O)C(O...  Rejected  Rejected  Rejected  ...   9  209.76    4            5
...                                                  ...       ...       ...       ...  ...   .     ...    .            .
20010  COC12C3NC3CN2C2=C(C1COC(=O)N)C(=O)C(=C(C2=O)C)...  Rejected  Rejected  Rejected  ...  10  375.28   24            7
20011              COC(=O)c1ccc(cc1F)N1CC(OC1=O)CNC(=S)C  Accepted  Accepted  Accepted  ...   1   67.87    6            3
20012                                 S=C=NC(=O)c1ccccc1  Accepted  Accepted  Accepted  ...   0   29.43    2            3
20013                                    S=C=NCCc1ccccc1  Accepted  Accepted  Accepted  ...   0   12.36    3            3
20014  CN1N=Cc2ccc(cc2)OP2(=NP3(=S)Oc4ccc(cc4)C=NN(C)...  Rejected  Rejected  Rejected  ...   0  223.36    6            6

>>> summary_2
              Rejected  Accepted
bRo5              3833     16182
Macro             5033     14982
Lipinski          4508     15507
Xu                4508     15507
NTD              13032      6983
Toxicophores     12380      7635
SureChemble       3013     17002

.. code-block:: python
	:linenos:

	import matplotlib.pyplot as plt
	
	f,ax = plt.subplots()

	labels = ['seven-filter','six-filter','five-filter','four-filter','three-filter','two-filter','one-filter','zero-filter']
	res = summary.Rejected_Num.value_counts().sort_index(ascending=False)
	sizes = res.values
	explode = (0,0,0,0,0,0,0,0.1)
	ax.pie(sizes,explode=explode,labels=labels,autopct='%1.1f%%',shadow=False,startangle=150)
	plt.show()  

.. figure:: /image/natural_ratio.png
	:width: 400px
	:align: center

From the above pie chart, we can see 88.6% of natural products has been detected at least one filter.

Next, we take The mol_0 (Oc1cc2[OH+][C-](c3ccccc3)C(=Cc2c(c1)O)O) as axample to show the radar plot

.. code-block:: python
	:linenos:
	
	from rdkit import Chem
	from scopy.visualize import pc_depict
		
	mol = Chem.MolFromSmiles('Oc1cc2[OH+][C-](c3ccccc3)C(=Cc2c(c1)O)O')
	prop_kws = {'MW':[0,500],'logP':[None,5],'nHD':[0,5],
	            'nHA':[0,10],'TPSA':[0,140]} #Linpiski's Rule and TPSA from Veber's Rule
	f = pc_depict.rule_radar(mol, prop_kws)
    

.. figure:: /image/natural_radar.png
	:width: 400px
	:align: center



Application 3. Screening Small Molecule Libarary
--------------------------------------------------
For small molecule libraries, in addition to screening for drug-like properties and toxicity, frequent hitters are also a focus of screening. Frequent hitters are compounds that show positive results in multiple HTS tests. Most frequently hit compounds are false positives through interference tests, while other promiscuous compounds that actually bind to most targets are difficult to make due to their poor selectivity.
We collected 408,598 molecules from Zelin's database. For drug-like screening, we use 3 drug-like rules (Lipinski, Pfizer and GSK) for filtering, and 3 toxicity endpoint filters (Potential_Electrophilic, Skin_Sensitization and LD50_oral Filter) are used for toxicity screening. Finally, 6 frequent hit compound filters (AlphaScreen_FHs, AlphaScreen_GST_FHs, AlphaScreen_HIS_FHs, Chelating, BMS and PAINS Filter) are used to filter potentially frequent hit compounds.

.. figure:: /image/fh_intro.png
	:width: 400px
	:align: center

.. code-block:: python
	:linenos:

	from rdkit import Chem
	import pandas as pd #This package should be installed
	from scopy.structure_alert.SmartsFilter import Filter
	from scopy.druglikeness.druglikeness import PC_properties, PC_rules

	data = pd.read_csv('Zelin.csv')
	data['mol'] = data.mol.map(lambda x: obsmitosmile(x))
	mols = [Chem.MolFromSmiles(smi) for smi in data.mol.values]

	rules = PC_rules(mols,n_jobs=-1,detail=True)
	screener = Filter(mols,n_jobs=20,detail=True)

	#=========================================================================
	# Drug-likeness rule
	#=========================================================================
	gsk = pd.DataFrame(rules.CheckGSKRule())
	pizer = pd.DataFrame(rules.CheckPfizerRule())
	lipiski = pd.DataFrame(rules.CheckLipinskiRule())

	summary_druglike = pd.DataFrame({'SMILES':data.mol.values,
	                                 'Lipinski':lipiski.Disposed,
	                                 'Pizer':pizer.Disposed,
	                                 'GSK':gsk.Disposed})

	summary_druglike['Rejected_Num'] = (summary_druglike=='Rejected').sum(axis=1)
	summary_2_druglike = (summary_druglike.iloc[:,1:]=='Rejected').sum(axis=0)
	summary_2_druglike.columns = ['Rejected']
	summary_2_druglike['Accepted'] = len(summary_druglike)-summary_2_druglike.Rejected.values

	#=========================================================================
	# Toxicity
	#=========================================================================
	ele = pd.DataFrame(screener.Check_Potential_Electrophilic())
	skin = pd.DataFrame(screener.Check_Skin_Sensitization())
	ld_50 = pd.DataFrame(screener.Check_LD50_Oral())

	summary_tox = pd.DataFrame({'SMILES':data.mol.values,
	                            'Potential_Electrophilic':ele.Disposed,
	                            'Skin_Sensitization':skin.Disposed,
	                            'LD50_Oral':ld_50.Disposed,
	                            })

	summary_tox['Rejected_Num'] = (summary_tox=='Rejected').sum(axis=1)
	summary_2_tox = (summary_tox.iloc[:,1:]=='Rejected').sum(axis=0)
	summary_2_tox.columns = ['Rejected']
	summary_2_tox['Accepted'] = len(summary_tox)-summary_2_tox.Rejected.values

	#=========================================================================
	# Frequent Hitters
	#=========================================================================
	alapha_fh = pd.DataFrame(screener.Check_AlphaScreen_FHs())
	gst = pd.DataFrame(screener.Check_AlphaScreen_GST_FHs())
	his = pd.DataFrame(screener.Check_AlphaScreen_HIS_FHs())
	che = pd.DataFrame(screener.Check_Chelating())
	bms = pd.DataFrame(screener.Check_BMS())
	pains = pd.DataFrame(screener.Check_PAINS())

	summary_fh = pd.DataFrame({'SMILES':data.mol.values,
	                           'AlphaScreen_FHs':alapha_fh.Disposed,
	                           'AlphaScreen_GST_FHs':gst.Disposed,
	                           'AlphaScreen_HIS_FHs':his.Disposed,
	                           'Chelating':che.Disposed,
	                           'BMS':bms.Disposed,
	                           'PAINS':bms.Disposed})

	summary_fh['Rejected_Num'] = (summary_fh=='Rejected').sum(axis=1)
	summary_2_fh = (summary_fh.iloc[:,1:]=='Rejected').sum(axis=0)
	summary_2_fh.columns = ['Rejected']
	summary_2_fh['Accepted'] = len(summary_fh)-summary_2_fh.Rejected.values

>>> summary_2_druglike
     Filter  Rejected  Accepted
0  Lipinski      8448    400149
1     Pizer     73193    335404
2       GSK    200192    208405
>>> 
>>> summary_2_tox
                    Filter  Rejected  Accepted
0  Potential_Electrophilic    111946    296651
1       Skin_Sensitization    310897     97700
2                LD50_Oral      9466    399131
>>> 
>>>summary_2_fh
                Filter  Rejected  Accepted
0      AlphaScreen_FHs         3    408594
1  AlphaScreen_GST_FHs      3500    405097
2  AlphaScreen_HIS_FHs      2954    405643
3            Chelating      6757    401840
4                  BMS     92553    316044
5                PAINS     19375    389222

.. figure:: /image/withdrawn_ratio.png
	:width: 400px
	:align: center

.. code-block:: python
	:linenos:

	import matplotlib.pyplot as plt

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
    

From the figure above, we can see that although only 37.6% of small molecules did not pass the drug-like rule filter, most of them were caused by the GSK rule. In the GSK rule, the molecular mass is limited to less than 400 in order to obtain better ADME/T properties, which is a bit strict for small molecules; nearly half of the molecules do not pass the toxicity filter, but some are made by "skin_sensitvity endpoint", which can be avoided by changing the route of administration; in the end, nearly three-quarters of the molecules passed the frequent hit compound filter.



Application 4. Screening Fragment
----------------------------------
Fragment-based drug discovery (FBDD) has become mainstream. FBDD has been widely applied in both academia and industry, as evidenced by the large number of papers from universities, non-profit research institutions, biotechnology companies and pharmaceutical companies.
We collected 2,877 small molecular fragments to build a library of molecular fragments. We use the drug-like rules Ro2 and Ro3 that apply to the fragments for filtering. We also use environmental toxicity filters for filtering.

.. code-block:: python
	:linenos:

	from rdkit import Chem
	import pandas as pd #This package should be installed
	from scopy.structure_alert.SmartsFilter import Filter
	from scopy.druglikeness.druglikeness import PC_rules

	data = pd.read_csv(file)
	mols = [Chem.MolFromSmiles(smi) for smi in data.mol.values]
	    
	rules = PC_rules(mols,n_jobs=24,detail=True)
	screener = Filter(mols,n_jobs=24,detail=True)

	#=========================================================================
	# Drug-likeness rule
	#=========================================================================
	ro2 = pd.DataFrame(rules.CheckRo2())
	ro3 = pd.DataFrame(rules.CheckRo3())

	summary_druglike = pd.DataFrame({'SMILES':data.mol.values,
	                                 'Ro2':ro2.Disposed,
	                                 'Ro3':ro3.Disposed,})

	summary_druglike['Rejected_Num'] = (summary_druglike=='Rejected').sum(axis=1)
	summary_2_druglike = pd.DataFrame((summary_druglike.iloc[:,1:-1]=='Rejected').sum(axis=0))
	summary_2_druglike.columns = ['Rejected']
	summary_2_druglike['Accepted'] = len(summary_druglike)-summary_2_druglike.Rejected.values

	#=========================================================================
	# Environment
	#=========================================================================
	nonbio = pd.DataFrame(screener.Check_NonBiodegradable())
	acuq = pd.DataFrame(screener.Check_Acute_Aquatic_Toxicity())

	summary_env = pd.DataFrame({'SMILES':data.mol.values,
	                        'NonBiodegradable':nonbio.Disposed,
	                        'Acute_Aquatic_Toxicity':acuq.Disposed})

	summary_env['Rejected_Num'] = (summary_env=='Rejected').sum(axis=1)
	summary_2_env = pd.DataFrame((summary_env.iloc[:,1:-1]=='Rejected').sum(axis=0))
	summary_2_env.columns = ['Rejected']
	summary_2_env['Accepted'] = len(summary_env)-summary_2_env.Rejected.values

>>> summary_druglike
  Filter  Rejected  Accepted
0    Ro2       833      2043
1    Ro3       188      2688
>>> 
>>> summary_env
                   Filter  Rejected  Accepted
0        NonBiodegradable      1147      1729
1  Acute_Aquatic_Toxicity      1450      1426

.. code-block:: python
	:linenos:

	import matplotlib.pyplot as plt

	f,axes = plt.subplots(1,2,figsize=(5*2,9))

	labels = ['two-filter','one-filter','zero-filter']

	res_druglike = summary_druglike.Rejected_Num.value_counts().sort_index(ascending=False)
	res_env = summary_env.Rejected_Num.value_counts().sort_index(ascending=False)

	axes[0].pie(res_druglike.values,explode=(0,0,0.1),
	            labels=labels,autopct='%1.1f%%',
	            shadow=False,startangle=150)
	axes[0].set_title('Drug-likeness Rule Filter')

	axes[1].pie(res_env.values,explode=(0,0,0.1),
	            labels=labels,autopct='%1.1f%%',
	            shadow=False,startangle=150)
	axes[1].set_title('Environmental Toxicity Filter')

	plt.show()   

.. figure:: /image/fragment_ratio.png
	:width: 400px
	:align: center