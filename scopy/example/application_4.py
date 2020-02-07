# -*- coding: utf-8 -*-

#Created on Thu Jan  2 14:30:43 2020
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    