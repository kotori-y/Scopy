# -*- coding: utf-8 -*-
"""
Created on Sat May 16 23:40:42 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


import os
from rdkit import Chem
from scopy.ScoConfig import DemoDir
from scopy.ScoTox import Toxfilter



def main(mols):
    tox = Toxfilter(mols, detail=False, showSMILES=True)
    
    tox.Check_Acute_Aquatic_Toxicity()
    tox.Check_Biodegradable()
    tox.Check_DNA_Binding()
    tox.Check_Developmental_Mitochondrial()
    tox.Check_Genotoxic_Carcinogenicity_Mutagenicity()
    tox.Check_Idiosyncratic()
    tox.Check_LD50_Oral()
    tox.Check_NTD()
    tox.Check_NonBiodegradable()
    tox.Check_NonGenotoxic_Carcinogenicity()
    tox.Check_Potential_Electrophilic()
    tox.Check_Reactive_Unstable_Toxic()
    tox.Check_Skin_Sensitization()
    tox.Check_SureChEMBL()
    tox.Check_Toxicophores()
    
    return 'Pass'
    
    
if '__main__' == __name__:
    file = os.path.join(DemoDir, '50.sdf')
    mols = Chem.SDMolSupplier(file)
    mols = [mol for mol in mols]
    print('====================== ScoTox ======================')
    res = main(mols)
    print(res)
    print('====================== ScoTox ======================')