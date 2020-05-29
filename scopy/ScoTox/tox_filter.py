# -*- coding: utf-8 -*-

#Created on Tue Jun 25 21:59:42 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me


from rdkit.Chem import AllChem as Chem
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
    
    
def Check_Acute_Aquatic_Toxicity(mol, detail=False, showSMILES=False):
    """Check molecule under Acute_Aquatic_Toxicity Filter,
    which presents a compound may cause toxicity to liquid(water).
    There are 99 SMARTS in this endpoint.
    
    Reference:
        (1) `Hermens, J. L. (1990)`_;
        (2) `Verhaar, Henk JM, Cees J. Van Leeuwen (1992)`_.
    
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Hermens, J. L. (1990):
        https://ehp.niehs.nih.gov/doi/abs/10.1289/ehp.9087219
    .. _Verhaar, Henk JM, Cees J. Van Leeuwen (1992):
        https://www.sciencedirect.com/science/article/pii/0045653592902805
        
    """
    Aquatic = Filter('Acute_Aquatic_Toxicity',detail, showSMILES)
    Aquatic.get_pattl()
    res = Aquatic.scan(mol)
    return res


def Check_Biodegradable(mol, detail=False, showSMILES=False):
    """
    Check molecule under Biodegradable Filter,
    which presents a compound may be Biodegradable.
    There are 9 SMARTS in this enpoint
    
    Reference:          
        (1) Existing Substances Program (CD-ROM), released April, 2004 (2003).
        
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    """
    Biodegradable = Filter('Biodegradable',detail, showSMILES)
    Biodegradable.get_pattl()
    res = Biodegradable.scan(mol)
    return res

def Check_Potential_Electrophilic(mol, detail=False, showSMILES=False):
    """
    Check molecule under Potential_Electrophilic Filter,
    which presents a compound would be more probably take part in electrophilic reaction, 
    and the electrophilic reaction is strongly assosiated with protein binding.
    There are 119 SMARTS in this endpoint.
    
    Reference:
        (1) `Enoch, S. J (2011)`_.
        
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    

    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Enoch, S. J (2011):
        https://www.tandfonline.com/doi/abs/10.3109/10408444.2011.598141
    
    """
    Potential_Electrophilic = Filter('Potential_Electrophilic',detail, showSMILES)
    Potential_Electrophilic.get_pattl()
    res = Potential_Electrophilic.scan(mol)
    return res


def Check_Genotoxic_Carcinogenicity_Mutagenicity(mol, detail=False, showSMILES=False):
    """
    Check molecule under Developmental_Mitochondrial Filter,
    which presents a compound may cause carcinogenicity or(and) mutagenicity through genotoxic mechanisms.
    There are 117 SMARTS in this endpoint.
       
    Reference:
        (1) `Benigni, Romualdo and Cecilia Bossa. (2008)`_;
        (2) `Ashby, John and Raymond W. Tennant. (1988)`_;
        (3) `Kazius, Jeroen, Ross McGuire, Roberta Bursi. (2005)`_;
        (4) `Bailey, Allan B. (2005)`_.
                
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Benigni, Romualdo and Cecilia Bossa. (2008):
        https://www.sciencedirect.com/science/article/pii/S1383574208000781
    .. _Ashby, John and Raymond W. Tennant. (1988):
        https://www.sciencedirect.com/science/article/pii/0165121888901140
    .. _Kazius, Jeroen, Ross McGuire, Roberta Bursi. (2005):
        https://pubs.acs.org/doi/abs/10.1021/jm040835a
    .. _Bailey, Allan B. (2005):
        https://www.sciencedirect.com/science/article/pii/S0273230005000553
    
    """
    Genotoxic_Carcinogenicity_Mutagenicity = Filter('Genotoxic_Carcinogenicity_Mutagenicity',detail, showSMILES)
    Genotoxic_Carcinogenicity_Mutagenicity.get_pattl()
    res = Genotoxic_Carcinogenicity_Mutagenicity.scan(mol)
    return res
    

def Check_LD50_Oral(mol, detail=False, showSMILES=False):
    """
    Check molecule under LD50_Oral Filter,
    which presents a compound may cause acute toxicity during oral administration;
    There are 20 SMARTS in this endpoint.
    
    Reference:
        (1) Tinkov OV (2019).
         
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    """
    LD50_Oral = Filter('LD50_Oral',detail, showSMILES)
    LD50_Oral.get_pattl()
    res = LD50_Oral.scan(mol)
    return res
    

def Check_NonBiodegradable(mol, detail=False, showSMILES=False):
    """
    Check molecule under NonBiodegradable Filter,
    which presents a compound may be non-biodegradable.
    There are 19 SMARTS in this enpoint.
    
    Reference:
        (1) Environment Canada. Existing Substances Program (CD-ROM), released April, 2004 (2003).
        
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    """
    NonBiodegradable = Filter('NonBiodegradable',detail, showSMILES)
    NonBiodegradable.get_pattl()
    res = NonBiodegradable.scan(mol)
    return res


def Check_NonGenotoxic_Carcinogenicity(mol, detail=False, showSMILES=False):
    """
    Check molecule under NonGenotoxic_Carcinogenicity Filter,
    which presents a compound may cause carcinogenicity or(and) mutagenicity through genotoxic mechanisms。
    There are 23 SMARTS in this endpoint.
    
    Reference:
        (1) `Benigni, Romualdo and Cecilia Bossa (2008)`_;
        (2) `Benigni, Romualdo, Cecilia Bossa and Olga Tcheremenskaia (2013)`_.
    
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Benigni, Romualdo and Cecilia Bossa (2008):
        https://www.sciencedirect.com/science/article/pii/S1383574208000781
    .. _Benigni, Romualdo, Cecilia Bossa and Olga Tcheremenskaia (2013):
        https://pubs.acs.org/doi/abs/10.1021/cr300206t
      
    """
    NonGenotoxic_Carcinogenicity = Filter('NonGenotoxic_Carcinogenicity',detail, showSMILES)
    NonGenotoxic_Carcinogenicity.get_pattl()
    res = NonGenotoxic_Carcinogenicity.scan(mol)
    return res


def Check_Skin_Sensitization(mol, detail=False, showSMILES=False):
    """
    Check molecule under Skin_Sensitization Filter,
    There are 155 SMARTS in this endpoint.
    
    Reference:
        (1) `Payne, M. P. and P. T. Walsh (1994)`_.
        (2) `Enoch, S. J., J. C. Madden and M. T. D. Cronin (2008)`_.
        (3) `Barratt, M. D (1994)`_.
            
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    

    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Payne, M. P. and P. T. Walsh (1994):
        https://pubs.acs.org/doi/pdf/10.1021/ci00017a019
    .. _Enoch, S. J., J. C. Madden and M. T. D. Cronin (2008):
        https://www.tandfonline.com/doi/abs/10.1080/10629360802348985
    .. _Barratt, M. D (1994):
        https://www.sciencedirect.com/science/article/pii/0887233394902445
    
    """
    Skin_Sensitization = Filter('Skin_Sensitization',detail, showSMILES)
    Skin_Sensitization.get_pattl()
    res = Skin_Sensitization.scan(mol)
    return res
    

def Check_SureChEMBL(mol, detail=False, showSMILES=False):
    """
    Check molecule under SureChEMBL Filter,
    which presents a compound would match one or more structural alerts and hence considered to have a MedChem unfriendly status.  
    There are 164 SMARTS in this endpoint.
    
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    

    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _SureChEMBL:
        https://www.surechembl.org/knowledgebase/169485
    
    """
    SureChEMBL = Filter('SureChEMBL',detail, showSMILES)
    SureChEMBL.get_pattl()
    res = SureChEMBL.scan(mol)
    return res
    

def Check_NTD(mol, detail=False, showSMILES=False):
    """
    Brenk has proposed 105 unwanted groups in HTS
    
    Reference:
        (1) `Brenk, Ruth (2008)`_.
          
    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    

    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
    
    .. _Brenk, Ruth (2008):
        https://onlinelibrary.wiley.com/doi/abs/10.1002/cmdc.200700139
        
    """        
    NTD = Filter('NTD',detail, showSMILES)
    NTD.get_pattl()
    res = NTD.scan(mol)
    return res


def Check_Toxicophores(mol, detail=False, showSMILES=False):
    """
    There 154 SMARTS in Toxicophres

    :param mol: The molecule to be scanned
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    

    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'Endpoint', otherwise 'MatchedAtoms' and 'MatchedNames' are also provided.
    :rtype: namedtuple
         
    """
    Toxicophores = Filter('Toxicophores',detail, showSMILES)
    Toxicophores.get_pattl()
    res = Toxicophores.scan(mol)
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
    res = Check_Acute_Aquatic_Toxicity(mol,detail=True)
    print(res)
    
    
    
    