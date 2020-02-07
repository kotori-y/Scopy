# -*- coding: utf-8 -*-

#Created on Tue Dec 31 13:12:54 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me
#
#♥I love Princess Zelda forever♥


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
    
    
import matplotlib.pyplot as plt
f,ax = plt.subplots()
   
labels = ['seven-filter','six-filter','five-filter','four-filter','three-filter','two-filter','one-filter','zero-filter']
res = summary.Rejected_Num.value_counts().sort_index(ascending=False)
sizes = res.values
explode = (0,0,0,0,0,0,0,0.1)
ax.pie(sizes,explode=explode,labels=labels,autopct='%1.1f%%',shadow=False,startangle=150)
plt.show()  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    