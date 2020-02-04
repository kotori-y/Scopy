# -*- coding: utf-8 -*-

#Created on Thu Dec 26 20:19:33 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me
#
#♥I love Princess Zelda forever♥


from rdkit import Chem
#from scopy.pretreat import pretreat
from scopy.druglikeness import druglikeness


class Druglikeness_example(object):
    """
    """
    def __init__(self,molfile):
        self.molfile = molfile
        suppl = Chem.SDMolSupplier(molfile)
        self.mols = [mol for mol in suppl if mol]

    def smallmol_druglikeness(self):
        """Here, we use two druglikeness to filter a database which suit for small molecules
        
        Two Rules:
            Lipinski Rule: MW<=500, logP<=5, nHD<=5, nHA <=10
            Egan Rule: 0<=tPSA<=132; -1<=logP<=6
        """
        print('start')
        pc = druglikeness.PC_rules(self.mols,n_jobs=-1)
        self.egan_res = pc.CheckEganRule()
        print('over')
        
        print('start')
        pc = druglikeness.PC_rules(self.mols,n_jobs=-1)
        self.lipinski_res = pc.CheckLipinskiRule()
        print('over')
        
    def macro_druglikeness(self):
        """
        """   
        print('start')
        pc = druglikeness.PC_rules(self.mols,n_jobs=-1)
        self.bRo5 = pc.CheckBeyondRo5()
        print('pver')
        
        print('start')
        pc = druglikeness.PC_rules(self.mols,n_jobs=-1)
        self.macro = pc.CheckOralMacrocycles()
        print('over')

    def building_block(self):
        """
        """
        print('start')
        pc = druglikeness.PC_rules(self.mols,n_jobs=-1)
        self.Ro2 = pc.CheckRo2()
        print('over')
        
        print('start')
        pc = druglikeness.PC_rules(self.mols,n_jobs=-1)
        self.Ro3 = pc.CheckRo3()
        print('over')


import time

start = time.process_time()
Filter = Druglikeness_example()
Filter.smallmol_druglikeness()
Filter.macro_druglikeness()
Filter.building_block()
end = time.process_time()
print(end - start)
    
    
    
    
    
    
    