# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:05:26 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

♥I love Princess Zelda forever♥
"""


from rdkit import Chem
import pandas as pd #This package should be installed
import matplotlib.pyplot as plt
from scopy.visualize import mcloud, highlight, pc_depict
from scopy.structure_alert.SmartsFilter import Filter
from scopy.druglikeness.druglikeness import PC_properties, PC_rules
import openbabel as ob
import os

def obsmitosmile(smiles):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "can")
    conv.SetOptions("K", conv.OUTOPTIONS)
    mol = ob.OBMol()
    conv.ReadString(mol, smiles)
    smile = conv.WriteString(mol)
    smile = smile.replace('\t\n', '')
    return smile


def getMol(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        pass
    else:
        mol = Chem.MolFromSmiles(obsmitosmile(smi))
    return mol

def main(file):
    """
    """
    def _getCol(dic):
        return list(dic.keys())
    #=========================================================================
    # Read library
    #=========================================================================
    data = pd.read_csv(file)
    mols = [getMol(smi) for smi in data.mol.values]
    name = os.path.splitext(file)[0]
       
    #=========================================================================
    # Instantiate, and 20 processers has been used.
    #=========================================================================
    # props = PC_properties(mols,n_jobs=20)
    rules = PC_rules(mols,n_jobs=4,detail=True)
    screener = Filter(mols,n_jobs=4,detail=True)
    
    # #=========================================================================
    # # Framework analyse
    # #=========================================================================
    # scount = mcloud.CountScaffold(mols,stype='Murcko')
    # scount = pd.DataFrame(scount,index=['Frequency']).T
    # #scount.to_csv('scount.txt',sep='\t',header=None)
    
    # #=========================================================================
    # # PC Properties.
    # #=========================================================================
    # prop_res = props.GetProperties(items=['MW','logP','nHA','nHD','TPSA']) #5 propperties are choossen
    # prop_res = pd.DataFrame(prop_res)
    
    #=========================================================================
    # Drug-likeness rule
    #=========================================================================
    print('===========DR===============')
    # gsk = pd.DataFrame(rules.CheckGSKRule())
    # pizer = pd.DataFrame(rules.CheckPfizerRule())
    lipinski = pd.DataFrame(rules.CheckLipinskiRule())
    bRo5 = pd.DataFrame(rules.CheckBeyondRo5())
       
    summary_druglike = pd.DataFrame({'SMILES': data.mol.values,
                                     'Lipinski': lipinski.Disposed,
                                     'bRo5': bRo5.Disposed})
    
    summary_druglike['Rejected_Num'] = (summary_druglike=='Rejected').sum(axis=1)
    summary_2_druglike = pd.DataFrame((summary_druglike.iloc[:,1:]=='Rejected').sum(axis=0))
    summary_2_druglike.columns = ['Rejected']
    summary_2_druglike['Accepted'] = len(summary_druglike)-summary_2_druglike.Rejected.values
    
    lipinski.to_csv(r'../res/{}/lipinski_{}.csv'.format(name,name),index=False)
    bRo5.to_csv(r'../res/{}/bRo5_{}.csv'.format(name,name),index=False)
    summary_druglike.to_csv(r'../res/{}/summary_druglike_{}.csv'.format(name,name),index=False)
    summary_2_druglike.to_csv(r'../res/{}/summary_2_druglike_{}.csv'.format(name,name),index=False)
    
    #=========================================================================
    # Toxicity
    #=========================================================================
    print('===========TOX===============')
    ele = pd.DataFrame(screener.Check_Potential_Electrophilic())
    skin = pd.DataFrame(screener.Check_Skin_Sensitization())
    ld_50 = pd.DataFrame(screener.Check_LD50_Oral())
    gene_res = pd.DataFrame(screener.Check_Genotoxic_Carcinogenicity_Mutagenicity()) #Genotoxic_Carcinogenicity_Mutagenicity
    nogene_res = pd.DataFrame(screener.Check_NonGenotoxic_Carcinogenicity()) #NonGenotoxic_Carcinogenicity
    
    ntd = pd.DataFrame(screener.Check_NTD())
    chembl = pd.DataFrame(screener.Check_SureChEMBL())
    toxico = pd.DataFrame(screener.Check_Toxicophores())
    
    summary_tox = pd.DataFrame({'SMILES':data.mol.values,
                                'Potential_Electrophilic':ele.Disposed,
                                'Skin_Sensitization':skin.Disposed,
                                'LD50_Oral':ld_50.Disposed,
                                'Genotoxic_Carcinogenicity_Mutagenicity':gene_res.Disposed,
                                'NonGenotoxic_Carcinogenicity':nogene_res.Disposed,
                                'NTD':ntd.Disposed,
                                'SureChEMBL':chembl.Disposed,
                                'Toxicophores':toxico.Disposed})
       
    summary_tox['Rejected_Num'] = (summary_tox=='Rejected').sum(axis=1)
    summary_2_tox = pd.DataFrame((summary_tox.iloc[:,1:]=='Rejected').sum(axis=0))
    summary_2_tox.columns = ['Rejected']
    summary_2_tox['Accepted'] = len(summary_tox)-summary_2_tox.Rejected.values
    
    ele.to_csv(r'../res/{}/ele_{}.csv'.format(name,name),index=False)
    skin.to_csv(r'../res/{}/skin_{}.csv'.format(name,name),index=False)
    ld_50.to_csv(r'../res/{}/LD50_{}.csv'.format(name,name),index=False)
    gene_res.to_csv(r'../res/{}/gene_{}.csv'.format(name,name),index=False)
    nogene_res.to_csv(r'../res/{}/nogene_{}.csv'.format(name,name),index=False)
    ntd.to_csv(r'../res/{}/NTD_{}.csv'.format(name,name),index=False)
    chembl.to_csv(r'../res/{}/SureChEMBL_{}.csv'.format(name,name),index=False)
    toxico.to_csv(r'../res/{}/toxicop_{}.csv'.format(name,name),index=False)
    
    summary_druglike.to_csv(r'../res/{}/summary_tox_{}.csv'.format(name,name),index=False)
    summary_2_druglike.to_csv(r'../res/{}/summary_2_tox_{}.csv'.format(name,name),index=False)
    
    #=========================================================================
    # Frequent Hitters
    #=========================================================================
    print('AlphaScreen_FHs...')
    alapha_fh = pd.DataFrame(screener.Check_AlphaScreen_FHs())
    
    print('AlphaScreen_GST_FHs...')
    gst = pd.DataFrame(screener.Check_AlphaScreen_GST_FHs())
    
    print('AlphaScreen_HIS_FHs...')
    his = pd.DataFrame(screener.Check_AlphaScreen_HIS_FHs())
    
    print('Chelating...')
    che = pd.DataFrame(screener.Check_Chelating())
    
    print('BMS...')
    bms = pd.DataFrame(screener.Check_BMS())
    
    print('PAINS...')
    pains = pd.DataFrame(screener.Check_PAINS())
    
    print('Luciferase_Inhibitory...')
    luc = pd.DataFrame(screener.Check_Luciferase_Inhibitory())
    
    print('Alarm_NMR')
    nmr = pd.DataFrame(screener.Check_Alarm_NMR())
    
    alapha_fh.to_csv(r'../res/{}/alapha_fh_{}.csv'.format(name,name),index=False)
    gst.to_csv(r'../res/{}/gst_{}.csv'.format(name,name),index=False)
    his.to_csv(r'../res/{}/his_{}.csv'.format(name,name),index=False)
    che.to_csv(r'../res/{}/chelat_fh_{}.csv'.format(name,name),index=False)
    bms.to_csv(r'../res/{}/bms_fh_{}.csv'.format(name,name),index=False)
    pains.to_csv(r'../res/{}/pains_fh_{}.csv'.format(name,name),index=False)
    luc.to_csv(r'../res/{}/lucifer_fh_{}.csv'.format(name,name),index=False)
    nmr.to_csv(r'../res/{}/nmr_fh_{}.csv'.format(name,name),index=False)
    
    summary_fh = pd.DataFrame({'SMILES':data.mol.values,
                               'AlphaScreen_FHs':alapha_fh.Disposed,
                               'Luciferase_Inhibitory':luc.Disposed,
                               'Alarm_NMR':nmr.Disposed,
                               'Chelating':che.Disposed,
                               'AlphaScreen_GST_FHs':gst.Disposed,
                               'AlphaScreen_HIS_FHs':his.Disposed,
                               'BMS':bms.Disposed,
                               'PAINS':pains.Disposed})
    summary_fh['Rejected_Num'] = (summary_fh=='Rejected').sum(axis=1)
    summary_2_fh = pd.DataFrame((summary_fh.iloc[:,1:]=='Rejected').sum(axis=0))
    summary_2_fh.columns = ['Rejected']
    summary_2_fh['Accepted'] = len(summary_fh)-summary_2_fh.Rejected.values
    
    summary_fh.to_csv(r'../res/{}/summary_fh_{}.csv'.format(name,name),index=False)
    summary_2_fh.to_csv(r'../res/{}/summary_2_fh_{}.csv'.format(name,name),index=False)
    
  





