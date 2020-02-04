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


def main():
    import os
    os.chdir(r'C:\Users\0720\Desktop\Project\scopy_ref\data\Scopy_databases\Molecules\Natural_product')
    
    data = pd.read_csv(r'Natural_product-NPS20015.csv')
    mols = [Chem.MolFromSmiles(smi) for smi in data.Canonical_SMILES.values]
    
    props = PC_properties(mols,n_jobs=-1)
    rules = PC_rules(mols,n_jobs=-1,detail=True)
    _filter = Filter(mols,n_jobs=4,detail=True)
    
    #=========================================================================
    # Drug-likeness rule
    #=========================================================================
    bRo5 = pd.DataFrame(rules.CheckBeyondRo5())
    Macro = pd.DataFrame(rules.CheckOralMacrocycles())
    lipiski = pd.DataFrame(rules.CheckLipinskiRule())
    xu = pd.DataFrame(rules.CheckLipinskiRule())
    
    bRo5.to_csv('bRo5.csv',index=False)
    Macro.to_csv('Macro.csv',index=False)
    lipiski.to_csv('lipiski.csv',index=False)
    xu.to_csv('xu.csv',index=False)
    
    #=========================================================================
    # Comprhenesive
    #=========================================================================
    ntd_res = pd.DataFrame(_filter.Check_NTD())
    toxci = pd.DataFrame(_filter.Check_Toxicophores())
    chemble = pd.DataFrame(_filter.Check_SureChEMBL())
    
    ntd_res.to_csv('NTD.csv',index=False)
    toxci.to_csv('Toxicophores.csv',index=False)
    chemble.to_csv('SureChEMBL.csv',index=False)
    
    #=========================================================================
    # Summary
    #=========================================================================
    
    out = pd.DataFrame({'mol':data.Canonical_SMILES.values,
                        'bRo5':bRo5.Disposed,
                        'Macro':Macro.Disposed,
                        'Lipinski':lipiski.Disposed,
                        'Xu':xu.Disposed,
                        'NTD':ntd_res.Disposed,
                        'Toxicophores':toxci.Disposed,
                        'SureChemble':chemble.Disposed,
                        'MW':props.CalculateMolWeight(),
                        'logP':props.CalculateLogP(),
                        'nHA':props.CalculateNumHAcceptors(),
                        'nHD':props.CalculateNumHDonors(),
                        'tPSA':props.CalculateTPSA(),
                        'nRot':props.CalculateNumRotatableBonds()})
    
    out.to_csv('Summary.csv',index=False)
    
    
    
if '__main__'==__name__:
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    