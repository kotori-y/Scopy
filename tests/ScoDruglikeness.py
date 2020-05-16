# -*- coding: utf-8 -*-
"""
Created on Sat May 16 22:04:51 2020

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
from scopy.ScoDruglikeness import PC_properties, PC_rules



def main(mols):
    props = PC_properties(mols)
    rules = PC_rules(mols, detail=True, showSMILES=True)
    
    return props.GetProperties(), rules.CheckLipinskiRule()
    
    
if '__main__' == __name__:
    file = os.path.join(DemoDir, '50.sdf')
    mols = Chem.SDMolSupplier(file)
    mols = [mol for mol in mols]
    print('====================== ScoDruglikeness ======================')
    props, rules = main(mols)
    print(props, rules)
    print('====================== ScoDruglikeness ======================')