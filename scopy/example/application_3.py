# -*- coding: utf-8 -*-

#Created on Tue Dec 31 14:08:35 2019
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    