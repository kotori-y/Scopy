# -*- coding: utf-8 -*-
"""
Created on Sun May 17 00:37:35 2020

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
from scopy.ScoRepresent import *



def main(mols):

    CalculateDaylight(mols)
    CalculateECFP(mols)
    CalculateEFG(mols)
    CalculateEState(mols)
    CalculateFCFP(mols)
    CalculateGhoseCrippen(mols)
    CalculateIFG(mols)
    CalculateMACCS(mols)
    CalculatePubChem(mols)
    CountCarbonScaffold(mols)
    CountMurckoFramework(mols)
    
    return 'Pass'

if '__main__' == __name__:    
    file = os.path.join(DemoDir, '50.sdf')
    mols = Chem.SDMolSupplier(file)
    mols = [mol for mol in mols]
    print('====================== ScoRepresent ======================')
    res = main(mols)
    print(res)
    print('====================== ScoRepresent ======================')
    
    
    
    
    
    
    
    
    