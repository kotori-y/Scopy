# -*- coding: utf-8 -*-
"""
Created on Fri May 15 21:47:43 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from multiprocessing import Pool
from rdkit import Chem

from ..structure_alert.FilterWithSmarts import Filter



class FHfilter(object):
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
        
    def Check_AlphaScreen_FHs(self):
        """
        Check molecule under Check_AlphaScreen_FHs Filter,
        which presents a compound may be alphascreen frequent hitters.
        There are 6 SMARTS in this endpoint.
            
        :return: a list of dictionary
        :rtype: list
        
        """
        AlphaScreen = Filter('AlphaScreen_FHs', self.detail, self.showSMILES)
        AlphaScreen.get_pattl()
        pool = Pool(self.n_jobs)
        AlphaScreen_FHs = pool.map_async(AlphaScreen.scan, self.mols).get()
        pool.close()
        pool.join()
        return AlphaScreen_FHs
        
    def Check_AlphaScreen_GST_FHs(self):
        """
        Check molecule under Check_AlphaScreen_GST_FHs Filter,
        which presents a compound may prevent GST/GSH interaction during HTS.
        There are 34 SMARTS in this endpoint.
        
        :return: a list of dictionary
        :rtype: list
        
        """
        GST = Filter('AlphaScreen_GST_FHs', self.detail, self.showSMILES)
        GST.get_pattl()
        pool = Pool(self.n_jobs)
        AlphaScreen_GST_FHs = pool.map_async(GST.scan, self.mols).get()
        pool.close()
        pool.join()
        return AlphaScreen_GST_FHs
    
    def Check_AlphaScreen_HIS_FHs(self):
        """
        Check molecule under Check_AlphaScreen_HIS_FHs Filter,
        which presents a compound prevents the binding of the protein His-tag moiety to nickel chelate.
        There are 19 SMARTS in this endpoint.
        
        :return: a list of dictionary
        :rtype: list
        
        """
        HIS = Filter('AlphaScreen_HIS_FHs', self.detail, self.showSMILES)
        HIS.get_pattl()
        pool = Pool(self.n_jobs)
        AlphaScreen_HIS_FHs = pool.map_async(HIS.scan, self.mols).get()
        pool.close()
        pool.join()
        return AlphaScreen_HIS_FHs
        
    def Check_Chelating(self):
        """
        Check molecule under Chelating Filter,
        which presents a compound may inhibit metalloproteins.
        Thers are 55 SMARTS in this endpoint
        
        :return: a list of dictionary
        :rtype: list
        
        """
        Che = Filter('Chelating', self.detail, self.showSMILES)
        Che.get_pattl()
        pool = Pool(self.n_jobs)
        Chelating = pool.map_async(Che.scan, self.mols).get()
        pool.close()
        pool.join()
        return Chelating
    
    def Check_Luciferase_Inhibitory(self):
        """
        There 3 SMARTS in Luciferase_Inhibitory Filter
        
        :return: a list of dictionary
        :rtype: list
        
        """
        luc = Filter('Luciferase_Inhibitory', self.detail, self.showSMILES)
        luc.get_pattl()
        pool = Pool(self.n_jobs)
        Luciferase_Inhibitory = pool.map_async(luc.scan, self.mols).get()
        pool.close()
        pool.join()
        return Luciferase_Inhibitory
                
    def Check_PAINS(self):
        """
        Check molecule under PAINS Filter,
        which presents a type of compounds tend to be hitted in HTS.
        There are 480 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        pai = Filter('Pains', self.detail, self.showSMILES)
        pai.get_pattl()
        pool = Pool(self.n_jobs)
        Pains = pool.map_async(pai.scan, self.mols).get()
        pool.close()
        pool.join()
        return Pains
               
    def Check_BMS(self):
        """
        Check molecule under BMS Filter.
        Pearce has proposed a Functional Group Compound Filters(FG Filters).
        The FG filters are consisted of two part, Exclusion FG filters and informational filters.
        Exclusion FG filters are those intended for compound removal from screening decks;
        Informational filters are useful for compound annotation.
        There are 176 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        bms = Filter('BMS', self.detail, self.showSMILES)
        bms.get_pattl()
        pool = Pool(self.n_jobs)
        BMS = pool.map_async(bms.scan, self.mols).get()
        pool.close()
        pool.join()
        return BMS
                
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
            
    def Check_Alarm_NMR(self):
        """
        There are 75 SMARTS in alarm_nmr
        
        :return: a list of dictionary
        :rtype: list
        
        """
        nmr = Filter('Alarm_NMR', self.detail, self.showSMILES)
        nmr.get_pattl()
        pool = Pool(self.n_jobs)
        Alarm_NMR = pool.map_async(nmr.scan, self.mols).get()
        pool.close()
        pool.join()
        return Alarm_NMR
            
    def Check_Frequent_Hitters(self):
        """
        There are 15 SMARTS in Frequent_Hitters

        :return: a list of dictionary
        :rtype: list
        
        """
        fh = Filter('Frequent_Hitters', self.detail, self.showSMILES)
        fh.get_pattl()
        pool = Pool(self.n_jobs)
        Frequent_Hitters = pool.map_async(fh.scan, self.mols).get()
        pool.close()
        pool.join()
        return Frequent_Hitters
        
    def Check_Aggregators(self):
        """
        There are 311 SMARTS in Aggregators
        
        :return: a list of dictionary
        :rtype: list
        
        """
        agg = Filter('Aggregators', self.detail, self.showSMILES)
        agg.get_pattl()
        pool = Pool(self.n_jobs)
        Aggregators = pool.map_async(agg.scan, self.mols).get()
        pool.close()
        pool.join()
        return Aggregators
            
    
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
    test = FHfilter(mols, n_jobs=4)
    res = test.Check_Luciferase_Inhibitory()
    print(res)