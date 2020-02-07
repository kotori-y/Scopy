# -*- coding: utf-8 -*-

#Created on Tue Dec 31 10:04:04 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me

#♥I love Princess Zelda forever♥


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



import matplotlib.pyplot as plt

f,axes = plt.subplots()

labels = ['six-filter','five-filter','four-filter','three-filter','two-filter','one-filter','zero-filter']
res = summary.Rejected_Num.value_counts().sort_index(ascending=False)
sizes = res.values
explode = (0,0,0,0,0,0,0.1)

ax.pie(sizes,explode=explode,labels=labels,autopct='%1.1f%%',shadow=False,startangle=150)
plt.show()  






