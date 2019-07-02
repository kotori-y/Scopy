# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:29:42 2019

@Author: CBDD Group, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com
@Blog: https://blog.moyule.me

"""

from scopy.StructureAlert import FliterWithSmarts
from scopy.StructureAlert import ComputeEFG
from rdkit.Chem import AllChem as Chem

def test_FliterWithSmarts():
    print('>>>>>>>>>>Testing FliterWithSmarts module\n\n')
    smis=['CCCCCC','CCC(C)CC','CC(C)CCC',
          'C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2','CCCCCN','c1ccccc1N']
    for idx,smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        print('{} The SMILES {}:'.format(str(idx+1),smi))
        res = FliterWithSmarts.Check_Acute_Aquatic_Toxicity(mol)
        print('\tThe all properties are: {}\n'.format(res))







if '__main__' == __name__:
    test_FliterWithSmarts()