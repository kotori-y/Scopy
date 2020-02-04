# -*- coding: utf-8 -*-

#Created on Tue Dec 31 10:04:04 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me

#♥I love Princess Zelda forever♥


from rdkit import Chem
import pandas as pd #This package should be installed
from scopy.structure_alert.SmartsFilter import Filter


def main():
    import os
    os.chdir(r'C:\Users\0720\Desktop\Project\scopy_ref\data\Scopy_databases\Molecules\Withdraw_drug')
    
    data = pd.read_csv(r"withdrawn_drug.csv")    
    mols = [Chem.MolFromSmiles(smi) for smi in data.mol.values]
    
    _filter = Filter(mols,n_jobs=4,detail=True)
    
    #=========================================================================
    # Toxicity
    #=========================================================================
    ele_res = pd.DataFrame(_filter.Check_Potential_Electrophilic())
    ld50_res = pd.DataFrame(_filter.Check_LD50_Oral())
    gene_res = pd.DataFrame(_filter.Check_Genotoxic_Carcinogenicity_Mutagenicity())
    no_gene = pd.DataFrame(_filter.Check_NonGenotoxic_Carcinogenicity())
    
    #=========================================================================
    # Comprhenesive
    #=========================================================================
    ntd_res = pd.DataFrame(_filter.Check_NTD())
    toxci = pd.DataFrame(_filter.Check_Toxicophores())
    chemble = pd.DataFrame(_filter.Check_SureChEMBL())
    
    ele_res.to_csv('Potential_Electrophilic.csv',index=False)
    ld50_res.to_csv('LD50_Oral.csv',index=False)
    gene_res.to_csv('Genotoxic_Carcinogenicity_Mutagenicity.csv',index=False)
    no_gene.to_csv('NonGenotoxic_Carcinogenicity.csv',index=False)
    
    ntd_res.to_csv('NTD.csv',index=False)
    toxci.to_csv('Toxicophores.csv',index=False)
    chemble.to_csv('SureChEMBL.csv',index=False)
    


if '__main__'==__name__:
    main()
    