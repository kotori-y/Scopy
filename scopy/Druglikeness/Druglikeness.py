# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 10:03:32 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

from scopy.Druglikeness import rulefilter
from scopy.Druglikeness import PropRadar
from rdkit import Chem


class OralSmallMol(object):
    """
    """  
    def __init__(self,mol):
        self.mol = mol
        
    def getresult(self):
        self.Egan = rulefilter.CheckEganRule(self.mol,detail=True)
        self.Veber = rulefilter.CheckVeberRule(self.mol,detail=True)
        self.Lipinski = rulefilter.CheckLipinskiRule(self.mol,detail=True)
#        oprea = rulefilter.CheckOpreaRule(mol)
        self.Xu = rulefilter.CheckXuRule(self.mol,detail=True)
#        CNS = rulefilter.CheckCNS(mol,detail=True)
#        kelder = rulefilter.CheckKelderRule(mol,detail=True)
#        reos = rulefilter.CheckREOS(mol,detail=True)
#        golden = rulefilter.CheckGoldenTriangle(mol,detail=True)
#        schneider = rulefilter.CheckSchneiderRule(mol,detail=True)            
        self.Result = {'egan':self.Egan,
                       'veber':self.Veber,
                       'lipinski':self.Lipinski,
                       'xu':self.Xu}
        
        
    def showXu(self):
        return PropRadar.VisualizeXu(XuRule=self.xu)
    
    def showLipinski(self):
        return PropRadar.VisualizeLipinski(LipinskiRule=self.lipinski)
    
#    def showEgan(self):
#        return PropRadar.VisualizeEgan(self.egan)
#


class OralBigMol(object):
    """
    """
    def __init__(self,mol):
        self.mol = mol
        
    def getresult(self):
        self.bRo5 = rulefilter.CheckBeyondRo5(self.mol,detail=True)
        self.OralMacro = rulefilter.CheckOralMacrocycles(self.mol,detail=True)
    
    def show(self):
        pass
    
    
class Toxicity(object):
    """
    """
    def __init__(self,mol):
        self.mol = mol
        
    def getresult(self):
        self.GSK = rulefilter.CheckGSKRule(self.mol,detail=True)
        self.Pfizer = rulefilter.CheckPfizerRule(self.mol,detail=True)
        self.Result = {'GSK':self.GSK,
                       'Pfizer':self.Pfizer,
                       }
        
    def showPfizer(self):
        return PropRadar.PfizerPositioning(self.Pfizer)
    

class Leadlikeness(object):
    pass


class Fragment(object):
    pass

    
    
    
    
if '__main__' == __name__:
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O',
            'CC(N)C(=O)[O-].[Na+]','CC(=O)OC1=CC=CC=C1C(=O)O']
    
    smi = smis[4]
    mol = Chem.MolFromSmiles(smi)
    demo = Toxicity(mol)
    demo.getresult()
    demo.showPfizer()
    
    
    
    
    
    
    
    
    
    
    
    