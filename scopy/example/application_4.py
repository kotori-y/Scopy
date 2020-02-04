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
import openbabel as ob
import pandas as pd #This package should be installed
from scopy.structure_alert.SmartsFilter import Filter
from scopy.druglikeness.druglikeness import PC_rules
import warnings
warnings.filterwarnings('ignore')
import datetime
import os
from multiprocessing import Pool
import re


def obsmitosmile(smi):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "can")
    conv.SetOptions("K", conv.OUTOPTIONS)
    mol = ob.OBMol()
    conv.ReadString(mol, smi)
    smile = conv.WriteString(mol)
    smile = smile.replace('\t\n', '')
    return smile

def getmol(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        return mol
    else:
        return Chem.MolFromSmiles(obsmitosmile(smi))

def main(folder):
    """
    """
    
    os.chdir(folder)
    
    print('+++++++++++++PREPARE++++++++++++++++')
    
    files = os.listdir()
    
    for file in files:
        
        data = pd.read_csv(file)
    #    data['mol'] = data.mol.map(lambda x: obsmitosmile(x))
        name = re.findall('\d\d\d\d',file)[0]
            
        ps = Pool(24)
        mols = ps.map_async(getmol, data.mol.values).get()
        ps.close()
        ps.join()
        
        rules = PC_rules(mols,n_jobs=24,detail=True)
        _filter = Filter(mols,n_jobs=24,detail=True)
        
        #=========================================================================
        # Drug-likeness rule
        #=========================================================================
        print('========== Drug-likeness rule ==========')
        
        ro2 = pd.DataFrame(rules.CheckRo2())
        print('>>> Ro2 finished')
        ro3 = pd.DataFrame(rules.CheckRo3())
        print('>>> Ro3 finished')
        
        ro2.insert(0,'SMILES',data.mol.values)
        ro3.insert(0,'SMILES',data.mol.values)
        
        ro2.to_csv('ro2_{}.csv'.format(name),index=False)
        ro3.to_csv('ro3_{}.csv'.format(name),index=False)
        
        summary = pd.DataFrame({'SMILES':data.mol.values,
                                'Ro2':ro2.Disposed,
                                'Ro3':ro3.Disposed,
                                })
        
        summary['Rejected_Num'] = (summary=='Rejected').sum(axis=1)
        summary_2 = pd.DataFrame((summary.iloc[:,1:-1]=='Rejected').sum(axis=0))
        summary_2.columns = ['Rejected']
        summary_2['Accepted'] = len(summary)-summary_2.Rejected.values
        
        summary.to_csv('{}_summary.csv'.format(name,name),index=False)
        summary_2.to_csv('{}_summary_2.csv'.format(name,name),index_label='Filter')
        
        
        #=========================================================================
        # Environment
        #=========================================================================
        print('\n\n========== Env ==========')
    
        nonbio = pd.DataFrame(_filter.Check_NonBiodegradable())
        print('>>> NonBiodegradable finished')
        acuq = pd.DataFrame(_filter.Check_Acute_Aquatic_Toxicity())
        print('>>> Acute_Aquatic_Toxicity finished')
        
        nonbio.insert(0,'SMILES',data.mol.values)
        acuq.insert(0,'SMILES',data.mol.values)
    
        nonbio.to_csv('nonBiode_{}.csv'.format(name),index=False)
        acuq .to_csv('aquatic_{}.csv'.format(name),index=False)
        
        summary = pd.DataFrame({'SMILES':data.mol.values,
                                'NonBiodegradable':nonbio.Disposed,
                                'Acute_Aquatic_Toxicity':acuq.Disposed,
                                })
        
        summary['Rejected_Num'] = (summary=='Rejected').sum(axis=1)
        summary_2 = pd.DataFrame((summary.iloc[:,1:-1]=='Rejected').sum(axis=0))
        summary_2.columns = ['Rejected']
        summary_2['Accepted'] = len(summary)-summary_2.Rejected.values
        
        summary.to_csv('{}_summary.csv'.format(name,name),index=False)
        summary_2.to_csv('{}_summary_2.csv'.format(name,name),index_label='Filter')
    
    
    
if '__main__'==__name__:
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    