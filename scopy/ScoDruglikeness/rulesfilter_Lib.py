# -*- coding: utf-8 -*-
"""
Created on Sat May 16 10:39:03 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from functools import partial
from multiprocessing import Pool
from rdkit import Chem
from . import rulesfilter, molproperty_Lib




def _GetSmi(mol):
    """
    Get the SMILES of molecule
    
    :param mols: molecule
    :type mols: rdkit.Chem.rdchem.Mol
    :return: The SMILES of molecule
    :rtype: string
    
    """
    return Chem.MolToSmiles(mol)


class PC_rules(object):
    """
    Here, we comdat the whole function that analyse PC-drived rules
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    
    """
    def __init__(self, mols, n_jobs=1, detail=False, showSMILES=False):
        self.mols = mols if type(mols) is not Chem.rdchem.Mol else [mols]
        self.detail = detail
        self.showSMILES = showSMILES
        self.n_jobs = n_jobs if n_jobs>=1 else None
        
    def CheckEganRule(self):
        """
        Bad or Good oral biovailability rule
    
        Reference:
            Egan, William J., Kenneth M. Merz, and John J. Baldwin. 
            J Med Chem, 43.21 (2000): 3867-3877.
        
        Rule details:
            0 <= TPSA <= 132
            -1 <= logP <=6
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        :return: the weight of molecule(contain hydrogen atoms)
        :rtype: list
        
        """
        fn = partial(rulesfilter.CheckEganRule, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
            
    def CheckVeberRule(self):
        """
        Bad or Good oral biovailability rule
        
        Reference:
            Veber, Daniel F., et al.
            Journal of medicinal chemistry 45.12 (2002): 2615-2623.
            
        Rule details:
            nRot <= 10
            TPSA <= 140
            nHB <= 12
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckVeberRule, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
        
    def CheckLipinskiRule(self):
        """
        Check molecular under Lipinski's rule
    
        Reference:
            Lipinski, Christopher A., et al.
            Advanced drug delivery reviews 23.1-3 (1997): 3-25.
        
        Rule details:
            MW <= 500
            logP <= 5
            nHD <= 5
            nHA <= 10
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckLipinskiRule, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
    
    def CheckBeyondRo5(self):
        """
        Check molecular under beyond Ro5
    
        Reference:
            Doak, Bradley C., et al.
            journal of medicinal chemistry 59.6 (2015): 2312-2327.
            
        Rule details:
            MW <= 1000
            -2 <= logP <= 10
            nHD <= 6
            nHA <= 15
            PSA <= 250
            nRot <= 20
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckBeyondRo5, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
    
    def CheckPfizerRule(self):
        """
        Check molecular under Rfizer Rule(3/75 Rule)
    
        Reference:
            Hughes, Jason D., et al. 
            Bioorganic & medicinal chemistry letters 18.17 (2008): 4872-4875.
            
        Rule details:
            logp > 3
            TPSA < 75
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckPfizerRule, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res        
    
    def CheckGSKRule(self):
        """
        Check molecular under GSK rule(4/400 Rule)
    
        Reference:
            Gleeson, M. Paul.
            Journal of medicinal chemistry 51.4 (2008): 817-834.
            
        Rule details:
            MW <= 400
            logP <= 4
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckGSKRule, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
        
    def CheckOralMacrocycles(self):
        """
        Check molecular under oral macrocycles rules
    
        Reference:
            Giordanetto, Fabrizio, and Jan Kihlberg.
            Journal of medicinal chemistry 57.2 (2013): 278-295.        
    
        Rule details:
            MW < 1000
            logP < 10
            nHD < 5
            TPSA < 250
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckOralMacrocycles, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
    
    def CheckOpreaRule(self):
        """
        Reference:
            Oprea, Tudor I.
            Journal of computer-aided molecular design 14.3 (2000): 251-264.
        
        Rules details:
            nRing >= 3
            nRig >= 18
            nRot >=6
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckOpreaRule, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
        
    def CheckGhoseRule(self):
        """
        Check molecular under Ghose rule
    
        Reference.:
            Ghose, Arup K., Vellarkad N. Viswanadhan, and John J. Wendoloski. 
            Journal of combinatorial chemistry 1.1 (1999): 55-68.
        
        Rules details:
            -0.4 < logP < 5.6
            160 < MW < 480
            40 < MR< 130
            20 < nAtom < 70
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckGhoseRule, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
        
    def CheckREOS(self):
        """
        Check molecular under REOS program
    
        Reference:
            Walters, W. Patrick, and Mark Namchuk.
            Nat Rev Drug Discov, 2.4 (2003): 259.
            
        Rule details:
            200 <= MW <= 500
            -5 <= logP <= 5
            nHD <= 5
            nHA <= 10
            nRot <= 8
            TPSA <= 150
            -4 <= fChar <= 4
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckREOS, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res
    
    def CheckGoldenTriangle(self):
        """
        Check molecular under 'Golden Triangle'
    
        Reference:
            Johnson, Ted W., Klaus R. Dress, and Martin Edwards.
            Bioorg Med Chem Lett, 19.19 (2009): 5560-5564.
            
        Rule details:
            200 <= MW <= 500
            -2 <= logD <= 5
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        pass
#        fn = partial(rulesfilter.CheckGoldenTriangle, detail=self.detail, showSMILES=self.showSMILES)
#        pool = Pool(self.n_jobs)
#        res = pool.map_async(fn, self.mols).get()
#        pool.close()
#        pool.join()
#        return res  
        
    def CheckXuRule(self):
        """
        Check molecular under Xu's rule
    
        Reference:
            
        
        Rule details:
            nhd <= 5
            nha <= 10
            3 <= rot <= 35
            1 <= nring <= 7
            10 <= nhev <= 50
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckXuRule, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res  
        
    def CheckSchneiderRule(self):
        """
        Check  molecular under Schneider rule
    
        Reference:
            Schneider, Nadine, et al.
            J Chem Inf Model, 48.3 (2008): 613-628.
            
        Rule details:
            mw > 230
            nhd > 0
            nha > 0
            nrot > 0
            nring > 0
            mr > 40
            functional groups > 0
            molecular volume > 191
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        pass
        
    def CheckRo4(self):
        """
        Referenece:
        
        Rule details:
            MW <= 400
            logP <= 4
            nHD <= 4
            nHA <= 8
            TPSA <= 120
        
        :param mols: molecules
        :type mols: Iterable
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: bool, optional
        
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckRo4, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res          
    
    def CheckRo3(self):
        """
        Check molecular under Ro3
    
        Ref.:
            Congreve, Miles, et al.
            Drug discovery today 19.8 (2003): 876-877.
            
        Rule details:
            MW <= 300
            -3 <= logP <= 3
            NHD <= 3
            NHA <= 6
            PSA <= 60
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckRo3, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res          
        
    def CheckRo2(self):
        """
        Check molecular under RO2
    
        Ref.:
            Goldberg, Frederick W., et al.
            Drug Discovery Today 20.1 (2015): 11-17.
            
        Rule details:
            MW <= 200
            Logp <= 2
            NHD <= 2
            NHA <= 4
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckRo2, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res 
    
    def CheckCNS(self):
        """
        Check mol under CNS
    
        Reference:
            Jeffrey, Phil, and Scott Summerfield.
            Neurobiol Dis, 37.1 (2010): 33-37.
            
        Rule details:
            135 <= weight <= 582
            -0.2 <= logP <= 6.1
            nha <= 5
            nhd <= 3
            3 <= TPSA <= 118
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckCNS, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res 
    
    def CheckRespiratory(self):
        """
        Check mol under Respiratory
    
        Reference:
            Ritchie, Timothy J., Christopher N. Luscombe, and Simon JF Macdonald. 
            J Chem Inf Model, 49.4 (2009): 1025-1032.
            
        Rule details:
            240<=MW<= 520
            -2.0<=logP<=4.7
            6<=nHB<=12
            51<=tPSA<=135
            3<=nRot<=8
            1<=nRing<=5
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `list`
        
        """
        fn = partial(rulesfilter.CheckRespiratory, detail=self.detail, showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
        return res 
    
    def CheckCustomizeRule(self, prop_kws, closed_interval=True):
        """
        You could customize the rule with mostly properties ypu want.
        
        :param prop_kws: the keys of dict are properties you want to check; the values should be a tuple or list with two elements,present the left- and right-bounded respectively.
        :type prop_kws: `dict`
        :param closed_interval: control whether using closed interval, defaults to True
        :type closed_interval: `bool`, optional
        :param detail: Control returning specific infomation or not, defaults to False
        :type detail: `bool`, optional
        :param showSMILES: Control returning SMILES or not, defaults to False
        :type showSMILES: bool, optional
        
        :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
        :rtype: `dict`
        """
        fn = partial(rulesfilter.CheckCustomizeRule, prop_kws=prop_kws, 
                     closed_interval=closed_interval, detail=self.detail, 
                     showSMILES=self.showSMILES)
        pool = Pool(self.n_jobs)
        res = pool.map_async(fn, self.mols).get()
        pool.close()
        pool.join()
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
            'S',
            'CCCCC(=O)[H]', #Biodegradable
            'C1=CN=C(C(=O)O)C=C1', #Chelating
            'C(OC)1=CC=C2OCC3OC4C=C(OC)C=CC=4C(=O)C3C2=C1',
            'C1=C2N=CC=NC2=C2N=CNC2=C1', #Genotoxic_Carcinogenicity_Mutagenicity
            'N(CC)(CCCCC)C(=S)N', #Idiosyncratic
            ]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    pc = PC_rules(mols,n_jobs=4,detail=True)
    res = pc.CheckLipinskiRule()
    print(len(res))