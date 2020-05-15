# -*- coding: utf-8 -*-
"""
Created on Thu May 14 16:54:59 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""

__doc__ == """This file is aim at comparing with 
FAF-Drugs4 (https://fafdrugs4.rpbs.univ-paris-diderot.fr/mobyle.html)
"""

from rdkit import Chem
from scopy.druglikeness.druglikeness import PC_properties, PC_rules


class FafCompare(PC_properties, PC_rules):
    
    def __init__(self, mols, n_jobs=1, detail=False, showSMILES=False):
        """
        """
        PC_properties.__init__(self, mols, n_jobs)
        PC_rules.__init__(self, mols, n_jobs, detail, showSMILES)
    
    def CalculateProperties(self):
        """
        """
        props = self.GetProperties()
        return props
    
    def ScreenRules(self):
        """
        """
        lipinski = self.CheckLipinskiRule()
        egan = self.CheckEganRule()
        veber = self.CheckVeberRule()
        pfizer = self.CheckPfizerRule()
        gsk = self.CheckGSKRule()
        oprea = self.CheckOpreaRule()
        ghose = self.CheckGhoseRule()
        xu = self.CheckXuRule()
        ro4 = self.CheckRo4()
        reos = self.CheckREOS()
        golden = self.CheckGoldenTriangle()
        
        items = ['lipinski', 'egan', 'veber',
                 'pfizer', 'gsk', 'oprea',
                 'ghose', 'xu', 'ro4',
                 'reos', 'golden']
        
        val = [lipinski, egan, veber,
               pfizer, gsk, oprea,
               ghose, xu, ro4,
               reos, golden]
        
        res = dict(zip(items, val))
        return res
    
    
if '__main__' == __name__:
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
            ]
    
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    test = FafCompare(mols, n_jobs=4)
    
    props = test.CalculateProperties()
    screen_res = test.ScreenRules()
