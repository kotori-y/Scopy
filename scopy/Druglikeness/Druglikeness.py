# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 10:03:32 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

from scopy.Druglikeness import CalculateProperty
from scopy.Druglikeness import CheckRule
from rdkit import Chem

class Likeness(object):
    
    def __init__(self,mol):
        self.mol = mol
        
    def Get_Properties(self):
        self.properties = CalculateProperty.GetProperties(self.mol)
    
    def Get_MW(self):
        self.MW = CalculateProperty.CalculateMolWeight(self.mol)
    
    def Get_nBond(self):
        self.nBond = CalculateProperty.CalculateNumBonds(self.mol)
           
    def Get_nAtom(self):
        self.nAtom = CalculateProperty.CalculateNumAtoms(self.mol)
    
    def Get_nHD(self):
        self.nHD = CalculateProperty.CalculateNumHDonors(self.mol)
        
    def Get_nHA(self):
        self.nHA = CalculateProperty.CalculateNumHAcceptors(self.mol)
    
    def Get_nHB(self):
        self.nHB = CalculateProperty.CalculateNumHyBond(self.mol)
    
    def Get_nHet(self):
        self.nHet = CalculateProperty.CalculateHeteroNumber(self.mol)
    
    def Get_nStero(self):
        
        
        
        
if '__main__' == __name__:
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O',
            'CC(N)C(=O)[O-].[Na+]','CC(=O)OC1=CC=CC=C1C(=O)O']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C',
          'CCCCCN','c1ccccc1N']
    smiring = ['C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2',
               'C1CCC2C3CC(C4CCCC5CCCCC45)CCC3CCC2C1',
               'C1CCC2(CCC3(CCCCC3)CC2)CC1']
    for smi in smis:
        mol = Chem.MolFromSmiles(smi)
        demo = Likeness(mol)
        demo.Get_nHD()
        print(demo.nHD)
    