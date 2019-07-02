# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 20:47:05 2019

@Author: CBDD Group, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com
@Blog: https://blog.moyule.me

"""

from scopy.Druglikeness import CalculateProperty
from scopy.Druglikeness import CheckRule
from rdkit.Chem import AllChem as Chem

def test_CalculateProperty():
    print('>>>>>>>>>>Testing CalculateProperty module\n\n')
    smis=['CCCCCC','CCC(C)CC','CC(C)CCC',
          'C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2','CCCCCN','c1ccccc1N']
    for idx,smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        print('{} The SMILES {}:'.format(str(idx+1),smi))
        res = CalculateProperty.GetProperties(mol)
        print('\tThe all properties are: {}\n'.format(res))
        
        

def test_CheckRule(detail=False):
    print('>>>>>>>>>>Testing CheckRule module\n\n')
    smis=['CCCCCC','CCC(C)CC','CC(C)CCC',
          'C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2','CCCCCN','c1ccccc1N']
    for idx,smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        print('{} The SMILES {}:'.format(str(idx+1),smi))
        res = CheckRule.CheckLipinskiRule(mol,detail=detail)
        print('\tThe all properties are: {}\n'.format(res))
        

if '__main__' == __name__:
    test_CheckRule()