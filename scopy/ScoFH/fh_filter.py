# -*- coding: utf-8 -*-

#Created on Tue Jun 25 21:59:42 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me


from rdkit import Chem
from .. import *



class Filter(object):
    """
    *Internal Use Only*
    the tool to check molecule(s) whether matched some unexpected endpoints.
    based on module SmartProcess
    
    :param endpoint: The name of endpoint for scanning molecule
    :type endpoint: str
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
     
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    """
    
    def __init__(self, endpoint, detail=False, showSMILES=False):
         self.endpoint = endpoint   
         self.detail = detail
         self.showSMILES = showSMILES
         
    def get_pattl(self):
        self.pattl = Loadpkl(self.endpoint)
    
    def scan(self,mol):
        return CheckWithSmarts(mol,
                               self.pattl, self.endpoint,
                               self.detail, self.showSMILES)
    

def Check_AlphaScreen_FHs(mol, detail=False, showSMILES=False):
    """
    Check molecule under Check_AlphaScreen_FHs Filter,
    which presents a compound may be alphascreen frequent hitters.
    There are 6 SMARTS in this endpoint.
    
    Reference:
        (1) `Schorpp, Kenji. (2014)`_.
        
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Schorpp, Kenji (2014):
        https://journals.sagepub.com/doi/full/10.1177/1087057113516861
        
    """
    AlphaScreen = Filter('AlphaScreen_FHs',detail, showSMILES)
    AlphaScreen.get_pattl()
    res = AlphaScreen.scan(mol)
    return res


def Check_AlphaScreen_GST_FHs(mol, detail=False, showSMILES=False):
    """
    Check molecule under Check_AlphaScreen_GST_FHs Filter,
    which presents a compound may prevent GST/GSH interaction during HTS.
    There are 34 SMARTS in this endpoint.
    
    References:
        (1) `Brenke, Jara K. (2016)`_.
    
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Brenke, Jara K. (2016):
        https://journals.sagepub.com/doi/abs/10.1177/1087057116639992
        
    """
    GST = Filter('AlphaScreen_GST_FHs',detail, showSMILES)
    GST.get_pattl()
    res = GST.scan(mol)
    return res


def Check_AlphaScreen_HIS_FHs(mol, detail=False, showSMILES=False):
    """
    Check molecule under Check_AlphaScreen_HIS_FHs Filter,
    which presents a compound prevents the binding of the protein His-tag moiety to nickel chelate.
    There are 19 SMARTS in this endpoint.
    
    Reference:
        (1) `Schorpp, Kenji. (2014)`_.
        
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
        
    .. _Schorpp, Kenji. (2014):
        https://journals.sagepub.com/doi/abs/10.1177/1087057113516861
        
    """
    HIS = Filter('AlphaScreen_HIS_FHs', detail, showSMILES)
    HIS.get_pattl()
    res = HIS.scan(mol)
    return res


def Check_Chelating(mol, detail=False, showSMILES=False):
    """
    Check molecule under Chelating Filter,
    which presents a compound may inhibit metalloproteins.
    Thers are 55 SMARTS in this endpoint
    
    Reference.:
        (1) `Agrawal, Arpita. (2010)`_.
    
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
        
    .. _Agrawal, Arpita. (2010):
        https://onlinelibrary.wiley.com/doi/abs/10.1002/cmdc.200900516
    
    """
    Chelating = Filter('Chelating', detail, showSMILES)
    Chelating.get_pattl()
    res = Chelating.scan(mol)
    return res


def Check_Luciferase_Inhibitory(mol, detail=False, showSMILES=False):
    """
    There 3 SMARTS in Luciferase_Inhibitory Filter
    
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    """
    Luciferase_Inhibitory = Filter('Luciferase_Inhibitory',detail, showSMILES)
    Luciferase_Inhibitory.get_pattl()
    res = Luciferase_Inhibitory.scan(mol)
    return res
    

def Check_PAINS(mol, detail=False, showSMILES=False):
    """
    Check molecule under PAINS Filter,
    which presents a type of compounds tend to be hitted in HTS.
    There are 480 SMARTS in this endpoint.
    
    Reference:
        (1) `Baell, Jonathan B. and Georgina A. Holloway (2010)`_.
    
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    

    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Baell, Jonathan B. and Georgina A. Holloway (2010):
        https://pubs.acs.org/doi/abs/10.1021/jm901137j
    
    """
    PAINS = Filter('Pains',detail, showSMILES)
    PAINS.get_pattl()
    res = PAINS.scan(mol)
    return res

    
def Check_BMS(mol, detail=False, showSMILES=False):
    """
    Check molecule under BMS Filter.
    Pearce has proposed a Functional Group Compound Filters(FG Filters).
    The FG filters are consisted of two part, Exclusion FG filters and informational filters.
    Exclusion FG filters are those intended for compound removal from screening decks;
    Informational filters are useful for compound annotation.
    There are 176 SMARTS in this endpoint.
    
    Reference:
        (1) `Pearce, Bradley C (2006)`_.
    
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    

    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Pearce, Bradley C (2006):
        https://pubs.acs.org/doi/abs/10.1021/ci050504m
    
    """
    BMS = Filter('BMS',detail, showSMILES)
    BMS.get_pattl()
    res = BMS.scan(mol)
    return res


def Check_Alarm_NMR(mol, detail=False, showSMILES=False):
    """
    There are 75 SMARTS in alarm_nmr 

    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    

    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
         
    """
    Alarm_NMR = Filter('Alarm_NMR',detail, showSMILES)
    Alarm_NMR.get_pattl()
    res = Alarm_NMR.scan(mol)
    return res
        

if __name__ == '__main__':
    smis = [
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
    mol = Chem.MolFromSmiles(smis[0])
    res = Check_PAINS(mol,detail=True)
    print(res)
    
    
    