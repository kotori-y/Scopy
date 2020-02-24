# -*- coding: utf-8 -*-

#Created on Mon Aug 12 14:48:23 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me


try:
	from .FilterWithSmarts import _Filter
except:
	from FilterWithSmarts import _Filter

from rdkit.Chem import AllChem as Chem
from multiprocessing import Pool

    
class Filter(object):
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
        self.mols = mols
        self.n_jobs = n_jobs
        self.detail = detail
        self.showSMILES = showSMILES
        
    def Check_Acute_Aquatic_Toxicity(self):
        """Check molecule under Acute_Aquatic_Toxicity Filter,
        which presents a compound may cause toxicity to liquid(water).
        There are 99 SMARTS in this endpoint.
        
        :return: a list of dictionary
        :rtype: list
        
        """
        Aquatic = _Filter('Acute_Aquatic_Toxicity', self.detail, self.showSMILES)
        Aquatic.get_pattl()
        pool = Pool(self.n_jobs)
        Acute_Aquatic_Toxicity = pool.map_async(Aquatic.scan, self.mols).get()
        pool.close()
        pool.join()
        return Acute_Aquatic_Toxicity
        
    def Check_AlphaScreen_FHs(self):
        """
        Check molecule under Check_AlphaScreen_FHs Filter,
        which presents a compound may be alphascreen frequent hitters.
        There are 6 SMARTS in this endpoint.
        
        :return: a list of dictionary
        :rtype: list
        
        """
        AlphaScreen = _Filter('AlphaScreen_FHs', self.detail, self.showSMILES)
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
        GST = _Filter('AlphaScreen_GST_FHs', self.detail, self.showSMILES)
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
        HIS = _Filter('AlphaScreen_HIS_FHs', self.detail, self.showSMILES)
        HIS.get_pattl()
        pool = Pool(self.n_jobs)
        AlphaScreen_HIS_FHs = pool.map_async(HIS.scan, self.mols).get()
        pool.close()
        pool.join()
        return AlphaScreen_HIS_FHs
        
    def Check_Biodegradable(self):
        """
        Check molecule under Biodegradable Filter,
        which presents a compound may be Biodegradable.
        There are 9 SMARTS in this enpoint
    
        :return: a list of dictionary
        :rtype: list
        
        """
        Biode = _Filter('Biodegradable', self.detail, self.showSMILES)
        Biode.get_pattl()
        pool = Pool(self.n_jobs)
        Biodegradable = pool.map_async(Biode.scan, self.mols).get()
        pool.close()
        pool.join()
        return Biodegradable
        
    def Check_Chelating(self):
        """
        Check molecule under Chelating Filter,
        which presents a compound may inhibit metalloproteins.
        Thers are 55 SMARTS in this endpoint
        
        :return: a list of dictionary
        :rtype: list
        
        """
        Che = _Filter('Chelating', self.detail, self.showSMILES)
        Che.get_pattl()
        pool = Pool(self.n_jobs)
        Chelating = pool.map_async(Che.scan, self.mols).get()
        pool.close()
        pool.join()
        return Chelating
        
    def Check_Developmental_Mitochondrial(self):
        """
        Check molecule under Developmental_Mitochondrial Filter,
        which presents a compound may casue Developmental Toxicity and Mitochondrial Toxicity.
        There are 12 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        Deve = _Filter('Developmental_Mitochondrial', self.detail, self.showSMILES)
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
        Geno = _Filter('Genotoxic_Carcinogenicity_Mutagenicity', self.detail, self.showSMILES)
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
        Idi = _Filter('Idiosyncratic', self.detail, self.showSMILES)
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
        ld = _Filter('LD50_Oral', self.detail, self.showSMILES)
        ld.get_pattl()
        pool = Pool(self.n_jobs)
        LD50_Oral = pool.map_async(ld.scan, self.mols).get()
        pool.close()
        pool.join()
        return LD50_Oral
                
    def Check_Luciferase_Inhibitory(self):
        """
        There 3 SMARTS in Luciferase_Inhibitory Filter
        
        :return: a list of dictionary
        :rtype: list
        
        """
        luc = _Filter('Luciferase_Inhibitory', self.detail, self.showSMILES)
        luc.get_pattl()
        pool = Pool(self.n_jobs)
        Luciferase_Inhibitory = pool.map_async(luc.scan, self.mols).get()
        pool.close()
        pool.join()
        return Luciferase_Inhibitory
                
    def Check_NonBiodegradable(self):
        """
        Check molecule under NonBiodegradable Filter,
        which presents a compound may be non-biodegradable.
        There are 19 SMARTS in this enpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        nonbi = _Filter('NonBiodegradable', self.detail, self.showSMILES)
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
        nonge = _Filter('NonGenotoxic_Carcinogenicity', self.detail, self.showSMILES)
        nonge.get_pattl()
        pool = Pool(self.n_jobs)
        NonGenotoxic_Carcinogenicity = pool.map_async(nonge.scan, self.mols).get()
        pool.close()
        pool.join()
        return NonGenotoxic_Carcinogenicity
                
    def Check_PAINS(self):
        """
        Check molecule under PAINS Filter,
        which presents a type of compounds tend to be hitted in HTS.
        There are 480 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        pai = _Filter('Pains', self.detail, self.showSMILES)
        pai.get_pattl()
        pool = Pool(self.n_jobs)
        Pains = pool.map_async(pai.scan, self.mols).get()
        pool.close()
        pool.join()
        return Pains
                
    def Check_Potential_Electrophilic(self):
        """
        Check molecule under Potential_Electrophilic Filter,
        which presents a compound would be more probably take part in electrophilic reaction, 
        and the electrophilic reaction is strongly assosiated with protein binding.
        There are 119 SMARTS in this endpoint.
    
        :return: a list of dictionary
        :rtype: list
        
        """
        ele = _Filter('Potential_Electrophilic', self.detail, self.showSMILES)
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
        reac = _Filter('Reactive_Unstable_Toxic', self.detail, self.showSMILES)
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
        skin = _Filter('Skin_Sensitization', self.detail, self.showSMILES)
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
        dna = _Filter('DNA_Binding', self.detail, self.showSMILES)
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
        chembl = _Filter('SureChEMBL', self.detail, self.showSMILES)
        chembl.get_pattl()
        pool = Pool(self.n_jobs)
        SureChEMBL = pool.map_async(chembl.scan, self.mols).get()
        pool.close()
        pool.join()
        return SureChEMBL
               
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
        bms = _Filter('BMS', self.detail, self.showSMILES)
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
        ntd = _Filter('NTD', self.detail, self.showSMILES)
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
        nmr = _Filter('Alarm_NMR', self.detail, self.showSMILES)
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
        fh = _Filter('Frequent_Hitters', self.detail, self.showSMILES)
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
        agg = _Filter('Aggregators', self.detail, self.showSMILES)
        agg.get_pattl()
        pool = Pool(self.n_jobs)
        Aggregators = pool.map_async(agg.scan, self.mols).get()
        pool.close()
        pool.join()
        return Aggregators
             
    def Check_Toxicophores(self):
        """
        There 154 SMARTS in Toxicophres
        
        :return: a list of dictionary
        :rtype: list
        
        """
        tox = _Filter('Toxicophores', self.detail, self.showSMILES)
        tox.get_pattl()
        pool = Pool(self.n_jobs)
        Toxicophores = pool.map_async(tox.scan, self.mols).get()
        pool.close()
        pool.join()
        return Toxicophores
        
        
if '__main__' == __name__:
    import time
    import warnings
    warnings.filterwarnings('ignore')

    
#    smis = [
#            'C1=CC=CC(C(Br)C)=C1',
#            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
#            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1',
#            'C1=NC(CCN)=CN1',
#            'C1CCCC(CCO)C1',
#            'C1=CC=C2N=C(O)C=CC2=C1',
#            'C(OC)1=C(C)C=C2OC[C@]([H])3OC4C(C)=C(OC)C=CC=4C(=O)[C@@]3([H])C2=C1C',
#            'C1=C2N=CC=NC2=C2N=CNC2=C1',
#            'C1=C(O)C=CC(O)=C1',
#            'CCC1(c2ccccc2)C(=O)NC(=O)NC1=O',
#            'N1=CN=CN=C1',
#            'C1=C2C=CC=CC2=CC2C=CC=CC1=2', #NonGenotoxic_Carcinogenicity
#            'C1=CC=C2C(=O)CC(=O)C2=C1', #Pains
#            'C1=CC=CC(COCO)=C1', #Potential_Electrophilic
#            'N1=NC=CN1C=O', #Promiscuity
#            'CC(=O)OC(=O)C1C=COC1', #Skin_Sensitization
#            'S',
#            'CCCCC(=O)[H]', #Biodegradable
#            'C1=CN=C(C(=O)O)C=C1', #Chelating
#            'C(OC)1=CC=C2OCC3OC4C=C(OC)C=CC=4C(=O)C3C2=C1',
#            'C1=C2N=CC=NC2=C2N=CNC2=C1', #Genotoxic_Carcinogenicity_Mutagenicity
#            'N(CC)(CCCCC)C(=S)N', #Idiosyncratic
#            ]
#    mols = [Chem.MolFromSmiles(smi) for smi in smis]*76
    from scopy import ScoConfig
    import os
    suppl = Chem.SDMolSupplier(os.path.join(ScoConfig.DemoDir,'760.sdf'))
    mols = (mol for mol in suppl if mol)
    print('------start------')
    start = time.process_time()
    F = Filter(mols, n_jobs=4, detail=True)
    res = F.Check_PAINS()
    end = time.process_time()
    print('Time cost: {}\nSamples: {}'.format((end-start), len(res)))
#    F.Check_Alarm_NMR()
#    print(F.Aggregators)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    