# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 10:03:32 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

try:
    from . import rulesfilter
except:
    import rulesfilter

class Druglikeness(object):
    """
    """
    def __init__(self, mols, detail=False, tolist=True):
        self.mols = mols
        self.detail = detail
        self.tolist = tolist
    
    def CheckEganRule(self):
        self.EganRule = map(lambda mol: rulesfilter.CheckEganRule(mol, detail=self.detail), self.mols)
        if self.tolist:
            self.EganRule = list(self.EganRule)
            
    def CheckVeberRule(self):
        self.VeberRule = list(map(lambda mol: rulesfilter.CheckEganRule(mol, detail=self.detail), self.mols))
        
    def CheckLipinskiRule(self):
        self.Lipinski = list(map(lambda mol: rulesfilter.CheckLipinskiRule(mol, detail=self.detail), self.mols))
    
    def CheckBeyondRo5(self):
        self.Lipinski = list(map(lambda mol: rulesfilter.CheckLipinskiRule(mol, detail=self.detail), self.mols))
    
    def CheckPfizerRule(self):
        self.PfizerRule = list(map(lambda mol: rulesfilter.CheckPfizerRule(mol, detail=self.detail), self.mols))
        
    def CheckGSKRule(self):
        self.GSKRule = list(map(lambda mol: rulesfilter.CheckGSKRule(mol, detail=self.detail), self.mols))
        
    def CheckOralMacrocycles(self):
        self.OralMacrocycles = list(map(lambda mol: rulesfilter.CheckOralMacrocycles(mol, detail=self.detail), self.mols))
        
    def CheckOpreaRule(self):
        self.OpreaRule = list(map(lambda mol: rulesfilter.CheckOpreaRule(mol, detail=self.detail), self.mols))
    
    def CheckGhoseRule(self):
        self.GhoseRule = list(map(lambda mol: rulesfilter.CheckGhoseRule(mol, detail=self.detail), self.mols))
     
    def CheckREOS(self):
        self.REOS = list(map(lambda mol: rulesfilter.CheckREOS(mol, self.detail), self.mols))
     
    def CheckGoldenTriangle(self):
        self.GoldenTriangle = list(map(lambda mol: rulesfilter.CheckGoldenTriangle(mol, self.detail), self.mols))
    
    def CheckXuRule(self):
        self.XuRule = list(map(lambda mol: rulesfilter.CheckXuRule(mol, self.detail), self.mols))
        
    def CheckSchneiderRule(self):
        self.SchneiderRule = list(map(lambda mol: rulesfilter.CheckSchneiderRule(mol, self.detail), mols))
        
    def CheckRo4(self):
        self.Ro4 = list(map(lambda mol: rulesfilter.CheckRo4(mol, self.detail), mols))
    
    def CheckRo3(self):
        self.Ro3 = list(map(lambda mol: rulesfilter.CheckRo3(mol, self.detail), mols))
        
    def CheckRo2(self):
        self.Ro2 = list(map(lambda mol: rulesfilter.CheckRo2(mol, self.detail), mols))
    
    def CheckCNS(self):
        self.CNS = list(map(lambda mol: rulesfilter.CheckCNS(mol, self.detail), mols))
    
    def CheckRespiratory(self):
        self.Respiratory = list(map(lambda mol: rulesfilter.CheckRespiratory(mol, self.detail), mols))
    
        
        
        
    
if '__main__' == __name__:
    from rdkit import Chem
    smis = [
            'C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1',
            'C1=NC(CCN)=CN1',
            'C1CCCC(CCO)C1',
            'C1=CC=C2N=C(O)C=CC2=C1',
            'C(OC)1=C(C)C=C2OC[C@]([H])3OC4C(C)=C(OC)C=CC=4C(=O)[C@@]3([H])C2=C1C',
            'C1=C2N=CC=NC2=C2N=CNC2=C1',
            'C1=C(O)C=CC(O)=C1',
            'CCC1(c2ccccc2)C(=O)NC(=O)NC1=O',
            'N1=CN=CN=C1',
            'C1=C2C=CC=CC2=CC2C=CC=CC1=2', #NonGenotoxic_Carcinogenicity
            'C1=CC=C2C(=O)CC(=O)C2=C1', #Pains
            'C1=CC=CC(COCO)=C1', #Potential_Electrophilic
            'N1=NC=CN1C=O', #Promiscuity
            'CC(=O)OC(=O)C1C=COC1', #Skin_Sensitization
            'S',
            'CCCCC(=O)[H]', #Biodegradable
            'C1=CN=C(C(=O)O)C=C1', #Chelating
            'C(OC)1=CC=C2OCC3OC4C=C(OC)C=CC=4C(=O)C3C2=C1',
            'C1=C2N=CC=NC2=C2N=CNC2=C1', #Genotoxic_Carcinogenicity_Mutagenicity
            'N(CC)(CCCCC)C(=S)N', #Idiosyncratic
            ]
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    rule = Druglikeness(mols, detail=True)
    rule.CheckRo3()
    res = rule.Ro3
    print(res)