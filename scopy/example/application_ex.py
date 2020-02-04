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
import openbabel as ob
import pandas as pd #This package should be installed
from scopy.structure_alert.SmartsFilter import Filter
#from scopy.druglikeness.druglikeness import PC_properties, PC_rules
import warnings
warnings.filterwarnings('ignore')
import datetime
import os
from multiprocessing import Pool


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

def main(folder,file):
    """
    """
    
    os.chdir(folder)
    
    print('+++++++++++++PREPARE++++++++++++++++')
    
    data = pd.read_csv(file)
#    data['mol'] = data.mol.map(lambda x: obsmitosmile(x))
    name,_ = os.path.splitext(file)
        
    ps = Pool(24)
    mols = ps.map_async(getmol, data.mol.values).get()
    ps.close()
    ps.join()
        
    with open(r'E:\student\kotori_y\example\Scopy_databases\out.log','a') as f_obj:
        f_obj.write(' '.join([datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), file, 'Start', str(len(mols))]))
    f_obj.close()
    
#    props = PC_properties(mols,n_jobs=-1)
#    rules = PC_rules(mols,n_jobs=24,detail=True)
    _filter = Filter(mols,n_jobs=24,detail=True)
    
    try:
        os.makedirs('{}_res/druglikeness'.format(name))
    except FileExistsError:
        pass
    try:
        os.makedirs('{}_res/fh'.format(name))
    except FileExistsError:
        pass
    try:
        os.makedirs('{}_res/tox'.format(name))
    except FileExistsError:
        pass
    
    
    print('+++++++++++++START++++++++++++++++')
    
#    #=========================================================================
#    # Tox& Comprehensive
#    #=========================================================================
    print('\n\n========== Toxicity Compounds ==========')

    ele = pd.DataFrame(_filter.Check_Potential_Electrophilic())
    print('>>> ele finished')
    ld_50 = pd.DataFrame(_filter.Check_LD50_Oral())
    print('>>> ld50 finished')
    gene = pd.DataFrame(_filter.Check_Genotoxic_Carcinogenicity_Mutagenicity())
    print('>>> Genotoxic finished')
    nogene = pd.DataFrame(_filter.Check_NonGenotoxic_Carcinogenicity())
    print('>>> NonGenotoxic finished')
    
    ntd = pd.DataFrame(_filter.Check_NTD())
    print('>>> NTD finished')
    toxci = pd.DataFrame(_filter.Check_Toxicophores())
    print('>>> Toxicophores finished')
    chemble = pd.DataFrame(_filter.Check_SureChEMBL())
    print('>>> sureChEMBL finished')
    
    ele.insert(0,'SMILES',data.mol.values)
    ld_50.insert(0,'SMILES',data.mol.values)
    gene.insert(0,'SMILES',data.mol.values)
    nogene.insert(0,'SMILES',data.mol.values)
    ntd.insert(0,'SMILES',data.mol.values)
    toxci.insert(0,'SMILES',data.mol.values)
    chemble.insert(0,'SMILES',data.mol.values)
    
    ele.to_csv('{}_res/tox/ele_{}.csv'.format(name,name),index=False)
    ld_50.to_csv('{}_res/tox/ld50_{}.csv'.format(name,name),index=False)
    ntd.to_csv('{}_res/tox/ntd_{}.csv'.format(name,name),index=False)
    toxci.to_csv('{}_res/tox/toxic_{}.csv'.format(name,name),index=False)
    chemble.to_csv('{}_res/tox/chembl_{}.csv'.format(name,name),index=False)
    gene.to_csv('{}_res/tox/_{}.csv'.format(name,name),index=False)
    nogene.to_csv('{}_res/tox/NonGenotoxic_{}.csv'.format(name,name),index=False)
    
    
    summary = pd.DataFrame({'SMILES':data.mol.values,
                            'Potential_Electrophilic':ele.Disposed,
                            'LD50_Oral':ld_50.Disposed,
                            'Genotoxic':gene.Disposed,
                            'NonGenotoxic':nogene.Disposed,
                            'NTD':ntd.Disposed,
                            'Toxicophores':toxci.Disposed,
                            'SureChEMBL':chemble.Disposed,
                            })
    
    summary['Rejected_Num'] = (summary=='Rejected').sum(axis=1)
    summary_2 = pd.DataFrame((summary.iloc[:,1:-1]=='Rejected').sum(axis=0))
    summary_2.columns = ['Rejected']
    summary_2['Accepted'] = len(summary)-summary_2.Rejected.values
    
    summary.to_csv('{}_res/tox/tox_summary_{}.csv'.format(name,name),index=False)
    summary_2.to_csv('{}_res/tox/tox_summary_2_{}.csv'.format(name,name),index_label='Filter')
    
    with open(r'E:\student\kotori_y\example\Scopy_databases\out.log','a') as f_obj:
        f_obj.write(' '.join([datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), file, 'Finished', str(len(mols))]))
    f_obj.close()
    
    
if '__main__'==__name__:
    mother = r'E:\student\kotori_y\example\Scopy_databases\Molecules'
    folders = ['Zelin']
    for folder in folders:
        folder = os.path.join(mother,folder)
        files = os.listdir(folder)
        for file in files:
            main(folder=folder, file=file)
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    