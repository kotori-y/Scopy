# -*- coding: utf-8 -*-
"""
Created on Fri May 15 22:27:53 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from multiprocessing import Pool
from rdkit import Chem


from ..ScoBase import Filter



class Toxfilter(object):
    """
    Here, we comdat the whole function that check endpoint retrieved from module FilterWithSmarts
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    """
    def __init__(self, mols, n_jobs=1, detail=False, showSMILES=False):
        """Initialization
        
        """
        self.mols = mols if type(mols) is not Chem.rdchem.Mol else [mols]
        self.n_jobs = n_jobs if n_jobs >=1 else None
        self.detail = detail
        self.showSMILES = showSMILES
        
    def Check_Acute_Aquatic_Toxicity(self):
        """Check molecule under Acute_Aquatic_Toxicity Filter,
        which presents a compound may cause toxicity to liquid(water).
        There are 99 SMARTS in this endpoint.
        
        :return: a list of dictionary
        :rtype: list
        
        """
        Aquatic = Filter('Acute_Aquatic_Toxicity', self.detail, self.showSMILES)
        Aquatic.get_pattl()
        pool = Pool(self.n_jobs)
        Acute_Aquatic_Toxicity = pool.map_async(Aquatic.scan, self.mols).get()
        pool.close()
        pool.join()
        return Acute_Aquatic_Toxicity
        
    def Check_Biodegradable(self):
        """
        Check molecule under Biodegradable Filter,
        which presents a compound may be Biodegradable.
        There are 9 SMARTS in this enpoint
    
        :return: a list of dictionary
        :rtype: list
        
        """
        Biode = Filter('Biodegradable', self.detail, self.showSMILES)
        Biode.get_pattl()
        pool = Pool(self.n_jobs)
        Biodegradable = pool.map_async(Biode.scan, self.mols).get()
        pool.close()
        pool.join()
        return Biodegradable
        
    def Check_Developmental_Mitochondrial(self):
        """
        Check molecule under Developmental_Mitochondrial Filter,
        which presents a compound may casue Developmental Toxicity and Mitochondrial Toxicity.
        There are 12 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        Deve = Filter('Developmental_Mitochondrial', self.detail, self.showSMILES)
        Deve.get_pattl()
        pool = Pool(self.n_jobs)
        Developmental_Mitochondrial = pool.map_async(Deve.scan, self.mols).get()
        pool.close()
        pool.join()
        return Developmental_Mitochondrial
                
    def Check_Genotoxic_Carcinogenicity_Mutagenicity(self):
        """
        Check molecule under Developmental_Mitochondrial Filter,
        which presents a compound may cause carcinogenicity or(and) mutagenicity through genotoxic mechanisms.
        There are 117 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        Geno = Filter('Genotoxic_Carcinogenicity_Mutagenicity', self.detail, self.showSMILES)
        Geno.get_pattl()
        pool = Pool(self.n_jobs)
        Genotoxic_Carcinogenicity_Mutagenicity = pool.map_async(Geno.scan, self.mols).get()
        pool.close()
        pool.join()
        return Genotoxic_Carcinogenicity_Mutagenicity
                
    def Check_Idiosyncratic(self):
        """
        Check molecule under Idiosyncratic Filter,
        which presents a compound may has diosyncratic toxicity.
        There are 35 SMARTS in this endpoint.
    
    
        :return: a list of dictionary
        :rtype: list
        
        """
        Idi = Filter('Idiosyncratic', self.detail, self.showSMILES)
        Idi.get_pattl()
        pool = Pool(self.n_jobs)
        Idiosyncratic = pool.map_async(Idi.scan, self.mols).get()
        pool.close()
        pool.join()
        return Idiosyncratic
        
    def Check_LD50_Oral(self):
        """
        Check molecule under LD50_Oral Filter,
        which presents a compound may cause acute toxicity during oral administration;
        There are 20 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        ld = Filter('LD50_oral', self.detail, self.showSMILES)
        ld.get_pattl()
        pool = Pool(self.n_jobs)
        LD50_Oral = pool.map_async(ld.scan, self.mols).get()
        pool.close()
        pool.join()
        return LD50_Oral
                
    def Check_NonBiodegradable(self):
        """
        Check molecule under NonBiodegradable Filter,
        which presents a compound may be non-biodegradable.
        There are 19 SMARTS in this enpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        nonbi = Filter('NonBiodegradable', self.detail, self.showSMILES)
        nonbi.get_pattl()
        pool = Pool(self.n_jobs)
        NonBiodegradable = pool.map_async(nonbi.scan, self.mols).get()
        pool.close()
        pool.join()
        return NonBiodegradable
             
    def Check_NonGenotoxic_Carcinogenicity(self):
        """
        Check molecule under NonGenotoxic_Carcinogenicity Filter,
        which presents a compound may cause carcinogenicity or(and) mutagenicity through genotoxic mechanisms。
        There are 23 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        nonge = Filter('NonGenotoxic_Carcinogenicity', self.detail, self.showSMILES)
        nonge.get_pattl()
        pool = Pool(self.n_jobs)
        NonGenotoxic_Carcinogenicity = pool.map_async(nonge.scan, self.mols).get()
        pool.close()
        pool.join()
        return NonGenotoxic_Carcinogenicity
                
    def Check_Potential_Electrophilic(self):
        """
        Check molecule under Potential_Electrophilic Filter,
        which presents a compound would be more probably take part in electrophilic reaction, 
        and the electrophilic reaction is strongly assosiated with protein binding.
        There are 119 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        ele = Filter('Potential_Electrophilic', self.detail, self.showSMILES)
        ele.get_pattl()
        pool = Pool(self.n_jobs)
        Potential_Electrophilic = pool.map_async(ele.scan, self.mols).get()
        pool.close()
        pool.join()
        return Potential_Electrophilic
        
    def Check_Reactive_Unstable_Toxic(self):
        """
        Check molecule under Reactive_Unstable_Toxic Filter.
        There are 335 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        reac = Filter('Reactive_Unstable_Toxic', self.detail, self.showSMILES)
        reac.get_pattl()
        pool = Pool(self.n_jobs)
        Reactive_Unstable_Toxic = pool.map_async(reac.scan, self.mols).get()
        pool.close()
        pool.join()
        return Reactive_Unstable_Toxic
            
    def Check_Skin_Sensitization(self):
        """
        Check molecule under Skin_Sensitization Filter,
        There are 155 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        skin = Filter('Skin_Sensitization', self.detail, self.showSMILES)
        skin.get_pattl()
        pool = Pool(self.n_jobs)
        Skin_Sensitization = pool.map_async(skin.scan, self.mols).get()
        pool.close()
        pool.join()
        return Skin_Sensitization
                
    def Check_DNA_Binding(self):
        """
        Check molecule under DNA_Binding Filter,
        There are 78 SMARTS in this endpoint.
        
        :return: a list of dictionary
        :rtype: list
        
        """
        dna = Filter('DNA_Binding', self.detail, self.showSMILES)
        dna.get_pattl()
        pool = Pool(self.n_jobs)
        DNA_Binding = pool.map_async(dna.scan, self.mols).get()
        pool.close()
        pool.join()
        return DNA_Binding
                
    def Check_SureChEMBL(self):
        """
        Check molecule under SureChEMBL Filter,
        which presents a compound would match one or more structural alerts and hence considered to have a MedChem unfriendly status.  
        There are 164 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        chembl = Filter('SureChEMBL', self.detail, self.showSMILES)
        chembl.get_pattl()
        pool = Pool(self.n_jobs)
        SureChEMBL = pool.map_async(chembl.scan, self.mols).get()
        pool.close()
        pool.join()
        return SureChEMBL
                               
    def Check_NTD(self, detail=False):
        """
        Brenk has proposed 105 unwanted groups in HTS
        
        :return: a list of dictionary
        :rtype: list
        
        """
        ntd = Filter('NTD', self.detail, self.showSMILES)
        ntd.get_pattl()
        pool = Pool(self.n_jobs)
        NTD = pool.map_async(ntd.scan, self.mols).get()
        pool.close()
        pool.join()
        return NTD
        
    def Check_Toxicophores(self):
        """
        There 154 SMARTS in Toxicophres
        
        :return: a list of dictionary
        :rtype: list
        
        """
        tox = Filter('Toxicophores', self.detail, self.showSMILES)
        tox.get_pattl()
        pool = Pool(self.n_jobs)
        Toxicophores = pool.map_async(tox.scan, self.mols).get()
        pool.close()
        pool.join()
        return Toxicophores
    

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
            'CCCCC(=O)[H]', #Biodegradable
            'C1=CN=C(C(=O)O)C=C1', #Chelating
            'C(OC)1=CC=C2OCC3OC4C=C(OC)C=CC=4C(=O)C3C2=C1',
            'C1=C2N=CC=NC2=C2N=CNC2=C1', #Genotoxic_Carcinogenicity_Mutagenicity
            'N(CC)(CCCCC)C(=S)N', #Idiosyncratic
            ]
    
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    test = Toxfilter(mols, n_jobs=4, detail=True)
    res = test.Check_SureChEMBL()
    print(res)