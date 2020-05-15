# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:27:57 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


import os

import pandas as pd
from rdkit import Chem

from scopy import ScoConfig
from scopy.druglikeness.druglikeness import PC_properties

sdf = os.path.join(ScoConfig.DemoDir, '760.sdf')
suppl = Chem.SDMolSupplier(sdf)

def CalculateProperties(mols):
    """
    """
    props = PC_properties(mols, n_jobs=4).GetProperties(items=['MW', 'logP', 'nHD', 'nHA', 'TPSA'])
    return props

def getmol(smi):
    mol = Chem.MolFromSmiles(smi)
    return mol

if '__main__'==__name__:
    data = pd.read_csv(r"C:\Users\0720\Desktop\Naturalproducts.csv")
    
    mols = [getmol(smi) for smi in data.SMILES]
    mols = [mol for mol in mols if mol]
    
    props = CalculateProperties(mols)