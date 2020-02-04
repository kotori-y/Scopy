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
from scopy.druglikeness.druglikeness import PC_properties, PC_rules


def obsmitosmile(smi):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "can")
    conv.SetOptions("K", conv.OUTOPTIONS)
    mol = ob.OBMol()
    conv.ReadString(mol, smi)
    smile = conv.WriteString(mol)
    smile = smile.replace('\t\n', '')
    return smile

def main(folder,file):
    """
    """
    import os
    os.chdir(folder)
    
    data = pd.read_csv(file)
    data['mol'] = data.mol.map(lambda x: obsmitosmile(x))
    name,_ = os.path.splitext(file)
    mols = [Chem.MolFromSmiles(smi) for smi in data.mol.values]
    
#    props = PC_properties(mols,n_jobs=-1)
    rules = PC_rules(mols,n_jobs=-1,detail=True)
    _filter = Filter(mols,n_jobs=20,detail=True)
    
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
    
    #=========================================================================
    # Drug-likeness rule
    #=========================================================================
    print('========== Drug-likeness rule ==========')
    
    gsk = pd.DataFrame(rules.CheckGSKRule())
    print('>>> gsk finished')
    pizer = pd.DataFrame(rules.CheckPfizerRule())
    print('>>> pizer finished')
    lipiski = pd.DataFrame(rules.CheckLipinskiRule())
    print('>>> lipinski finished')
    
    gsk.insert(0,'SMILES',data.mol.values)
    pizer.insert(0,'SMILES',data.mol.values)
    lipiski.insert(0,'SMILES',data.mol.values)
    
    gsk.to_csv('{}_res/druglikeness/gsk_{}.csv'.format(name,name),index=False)
    lipiski.to_csv('{}_res/druglikeness/lipiski_{}.csv'.format(name,name),index=False)
    pizer.to_csv('{}_res/druglikeness/pizer_{}.csv'.format(name,name),index=False)
    
    summary = pd.DataFrame({'SMILES':data.mol.values,
                            'Lipinski':lipiski.Disposed,
                            'Pizer':pizer.Disposed,
                            'GSK':gsk.Disposed})
    
    summary['Rejected_Num'] = (summary=='Rejected').sum(axis=1)
    summary_2 = (summary.iloc[:,1:]=='Rejected').sum(axis=0)
    summary_2.columns = ['Rejected']
    summary_2['Accepted'] = len(summary)-summary_2.Rejected.values
    
    summary.to_csv('{}_res/druglikeness/druglikeness_summary_{}.csv',index=False)
    summary_2.to_csv('{}_res/druglikeness/druglikeness_summary_2_{}.csv',index_label='Filter')
    
    #=========================================================================
    # Frequent Hitters
    #=========================================================================
    print('\n\n========== Frequent Hitters ==========')
    
    alapha_fh = pd.DataFrame(_filter.Check_AlphaScreen_FHs())
    print('>>> alapha_fh finished')
    gst = pd.DataFrame(_filter.Check_AlphaScreen_GST_FHs())
    print('>>> gst finished')
    his = pd.DataFrame(_filter.Check_AlphaScreen_HIS_FHs())
    print('>>> his finished')
    che = pd.DataFrame(_filter.Check_Chelating())
    print('>>> che finished')
    bms = pd.DataFrame(_filter.Check_BMS())
    print('>>> bms finished')
    pains = pd.DataFrame(_filter.Check_PAINS())
    print('>>> PAINS finished')
    
    alapha_fh.insert(0,'SMILES',data.mol.values)
    gst.insert(0,'SMILES',data.mol.values)
    his.insert(0,'SMILES',data.mol.values)
    che.insert(0,'SMILES',data.mol.values)
    bms.insert(0,'SMILES',data.mol.values)
    pains.insert(0,'SMILES',data.mol.values)
    
    alapha_fh.to_csv('{}_res/fh/alapha_fh_{}.csv'.format(name,name),index=False)
    gst.to_csv('{}_res/fh/gst_{}.csv'.format(name,name),index=False)
    his.to_csv('{}_res/fh/his_{}.csv'.format(name,name),index=False)
    che.to_csv('{}_res/fh/che_{}.csv'.format(name,name),index=False)
    bms.to_csv('{}_res/fh/bms_{}.csv'.format(name,name),index=False)
    pains.to_csv('{}_res/fh/pains_{}.csv'.format(name,name),index=False)
    
    summary = pd.DataFrame({'SMILES':data.mol.values,
                            'AlphaScreen_FHs':alapha_fh.Disposed,
                            'AlphaScreen_GST_FHs':gst.Disposed,
                            'AlphaScreen_HIS_FHs':his.Disposed,
                            'Chelating':che.Disposed,
                            'BMS':bms.Disposed,
                            'PAINS':bms.Disposed})
    
    summary['Rejected_Num'] = (summary=='Rejected').sum(axis=1)
    summary_2 = (summary.iloc[:,1:]=='Rejected').sum(axis=0)
    summary_2.columns = ['Rejected']
    summary_2['Accepted'] = len(summary)-summary_2.Rejected.values
    
    summary.to_csv('{}_res/fh/fh_summary_{}.csv',index=False)
    summary_2.to_csv('{}_res/fh/fh_summary_2_{}.csv',index_label='Filter')
    
    #=========================================================================
    # Tox
    #=========================================================================
    print('\n\n========== Toxicity Compounds ==========')

    ele = pd.DataFrame(_filter.Check_Potential_Electrophilic())
    print('>>> ele finished')
    skin = pd.DataFrame(_filter.Check_Skin_Sensitization())
    print('>>> gst finished')
    ld_50 = pd.DataFrame(_filter.Check_LD50_Oral())
    print('>>> ld50 finished')
    
    ele.insert(0,'SMILES',data.mol.values)
    skin.insert(0,'SMILES',data.mol.values)
    ld_50.insert(0,'SMILES',data.mol.values)

    ele.to_csv('{}_res/tox/ele_{}.csv'.format(name,name),index=False)
    skin.to_csv('{}_res/tox/skin_{}.csv'.format(name,name),index=False)
    ld_50.to_csv('{}_res/tox/ld50_{}.csv'.format(name,name),index=False)
    
    summary = pd.DataFrame({'SMILES':data.mol.values,
                            'Potential_Electrophilic':ele.Disposed,
                            'Skin_Sensitization':skin.Disposed,
                            'LD50_Oral':ld_50.Disposed,
                            })
    
    summary['Rejected_Num'] = (summary=='Rejected').sum(axis=1)
    summary_2 = (summary.iloc[:,1:]=='Rejected').sum(axis=0)
    summary_2.columns = ['Rejected']
    summary_2['Accepted'] = len(summary)-summary_2.Rejected.values
    
    summary.to_csv('{}_res/tox/tox_summary_{}.csv',index=False)
    summary_2.to_csv('{}_res/tox/tox_summary_2_{}.csv',index_label='Filter')
    
    
    
if '__main__'==__name__:
    main()
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    