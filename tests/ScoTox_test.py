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
    
    return tox.Check_LD50_Oral()
    
    
if '__main__' == __name__:
    file = os.path.join(DemoDir, '50.sdf')
    mols = Chem.SDMolSupplier(file)
    mols = [mol for mol in mols]
    print('====================== ScoTox ======================')
    res = main(mols)
    print(res)
    print('====================== ScoTox ======================')