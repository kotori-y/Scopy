# -*- coding: utf-8 -*-

#Created on Thu Dec 26 22:22:12 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me
#
#♥I love Princess Zelda forever♥


from rdkit import Chem
from scopy.structure_alert import SmartsFilter

class FrquentHitters_example(object):
    """
    """
    def __init__(self,molfile):
        molfile = molfile
        self.suppl = Chem.SDMolSupplier(molfile)

    def Check_PAINS(self):
        """
        """
        mols = (mol for mol in self.suppl if mol)
        self.Filter = SmartsFilter.Filter(mols, n_jobs=20)
        print('start')
        self.Pains_res = self.Filter.Check_PAINS()
        print('over')
        
    def Check_BMS(self):
        """
        """
        mols = (mol for mol in self.suppl if mol)
        self.Filter = SmartsFilter.Filter(mols, n_jobs=20)
        print('start')
        self.BMS_res = self.Filter.Check_BMS()
        print('over')
        
    def Check_Chelating(self):
        """
        """
        mols = (mol for mol in self.suppl if mol)
        self.Filter = SmartsFilter.Filter(mols, n_jobs=20)
        print('start')
        self.Chelating_res = self.Filter.Check_Chelating()
        print('over')
    
    def Check_Toxicity(self):
        """
        """
        mols = (mol for mol in self.suppl if mol)
        self.Filter = SmartsFilter.Filter(mols, n_jobs=20)
        print('start')
        self.Toxicity_res = self.Filter.Check_Toxicophores()
        print('over')
    
    def Check_Genotoxic(self):
        """
        """
        mols = (mol for mol in self.suppl if mol)
        self.Filter = SmartsFilter.Filter(mols, n_jobs=20)
        print('start')
        self.Genotoxic_res = self.Filter.Check_Genotoxic_Carcinogenicity_Mutagenicity()
        print('over')
        
    def Check_NonGenotoxic(self):
        """
        """
        mols = (mol for mol in self.suppl if mol)
        self.Filter = SmartsFilter.Filter(mols, n_jobs=20)
        print('start')
        self.NonGenotoxic_res = self.Filter.Check_NonGenotoxic_Carcinogenicity()
        print('over')
        
    def Check_NonBiodegradable(self):
        """
        """
        mols = (mol for mol in self.suppl if mol)
        self.Filter = SmartsFilter.Filter(mols, n_jobs=20)
        print('start')
        self.NonBiodegradable_res = self.Filter.Check_NonBiodegradable()
        print('over')
    
    def Check_Acuatic(self):
        """
        """
        mols = (mol for mol in self.suppl if mol)
        self.Filter = SmartsFilter.Filter(mols, n_jobs=20)
        print('start')
        self.Aquatic_res = self.Filter.Check_Acute_Aquatic_Toxicity()
        print('over')

import time

start = time.process_time()
F = FrquentHitters_example(r"C:\Anaconda\envs\my-rdkit-env\Lib\site-packages\scopy\data\Demo\760.sdf")
F.Check_BMS()
F.Check_Chelating()
F.Check_Toxicity()
F.Check_Genotoxic()
F.Check_NonBiodegradable()
F.Check_NonGenotoxic()
F.Check_Acuatic()

end = time.process_time()
print(end-start)