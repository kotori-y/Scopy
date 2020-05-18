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
    try:
        CalculateDaylight(mols)
    except:
        raise 
    else:
        print('>>> Daylight Pass')
        
    try:
        CalculateECFP(mols)
    except:
        raise 
    else:
        print('>>> ECFP Pass')
        
    try:
        CalculateEFG(mols)
    except:
        raise 
    else:
        print('>>> EFG Pass')
        
    try:
        CalculateEState(mols)
    except :
        raise 
    else:
        print('>>> EState Pass')
        
    try:
        CalculateFCFP(mols)
    except :
        raise 
    else:
        print('>>> FCFP Pass')
        
    try:
        CalculateGhoseCrippen(mols)
    except :
        raise 
    else:
        print('>>> GC Pass')
        
    try:
        CalculateIFG(mols)
    except :
        raise 
    else:
        print('>>> IFG Pass')
        
        
    try:
        CalculateMACCS(mols)
    except :
        raise 
    else:
        print('>>> MACCS Pass')
    
    try:
        CalculatePubChem(mols)
    except :
        raise 
    else:
        print('>>> PubChem Pass')
        
    try:
        CountCarbonScaffold(mols)
    except :
        raise 
    else:
        print('>>> Scaffold Pass')
        
    try:
        CountMurckoFramework(mols)
    except :
        raise 
    else:
        print('>>> FrameWork Pass')



if '__main__' == __name__:    
    file = os.path.join(DemoDir, '50.sdf')
    mols = Chem.SDMolSupplier(file)
    mols = [mol for mol in mols]
    print('====================== ScoRepresent ======================')
    main(mols)
    #print(res)
    print('====================== ScoRepresent ======================')
    
    
    
    
    
    
    
    
    