# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 21:59:42 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""
from __future__ import print_function
from rdkit.Chem import AllChem as Chem
import SmartProcess


class _Filter(object):
    """
    *Internal Use Only*
    the tool to check molecule(s) whether matched some unexpected endpoints.
    based on module SmartProcess
    
    Parameters:
    -----------
    endpoint: string
        The name of endpoint for scanning molecule(s).
    detail: bool
        Whether return more specific infomation.
        
    Return:
    -----------
    A namedtuple, 
    """
    
    def __init__(self,endpoint,detail):
         self.endpoint = endpoint   
         self.detail = detail
    
    def get_pattl(self):
        self.pattl = SmartProcess._Loadpkl(self.endpoint)
    
    def scan(self,mol):
        return SmartProcess._CheckWithSmarts(mol,self.pattl,self.endpoint,self.detail)
    
    
def Check_Acute_Aquatic_Toxicity(mols, detail=False):
    """
    Ref.:
    -----------
    (1) Hermens, J. L.
    Environ Health Perspect, 87 (1990): 219-225.
    (2) Verhaar, Henk JM, Cees J. Van Leeuwen, and Joop LM Hermens.
    Chemosphere, 25.4 (1992): 471-491.
    
    Brief:
    ----------- 
    The endpoint, Acute_Aquatic_Toxicity presents
    a compound may cause toxicity to liquid(water).
    There are 99 SMARTS in this endpoint.
   
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a list of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Acute_Aquatic_Toxicity(mols, detail=True)
    """
    Aquatic = _Filter('Acute_Aquatic_Toxicity',detail)
    Aquatic.get_pattl()
    res = list(map(Aquatic.scan,mols))
    return res


def Check_AlphaScreen_FHs(mols, detail=False):
    """
    Ref.:
    -----------
    Schorpp, Kenji, et al.
    J Biomol Screen, 19.5 (2014): 715-726.
        
    Brief:
    -----------
    The endpoint, AlphaScreen_FHs presents
    a compound may be alphascreen frequent hitters
    There are 6 SMARTS in this endpoint
   
    Parameters:
    -----------         
    mol: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_AlphaScreen_FHs(mols, detail=True)
    """
    AlphaScreen = _Filter('AlphaScreen_FHs',detail)
    AlphaScreen.get_pattl()
    res = list(map(AlphaScreen.scan,mols))
    return res


def Check_AlphaScreen_GST_FHs(mols, detail=False):
    """
    Ref.:
    -----------
    Brenke, Jara K., et al.
    J Biomol Screen, 21.6 (2016): 596-607.
    
    Brief:
    -----------
    The endpoint GST_FHs presents
    a compound may prevent GST/GSH interaction during HTS.
    There are 34 SMARTS in this endpoint
        
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_AlphaScreen_GST_FHs(mols, detail=True)
    """
    GST = _Filter('AlphaScreen_GST_FHs',detail)
    GST.get_pattl()
    res = list(map(GST.scan,mols))
    return res


def Check_AlphaScreen_HIS_FHs(mols, detail=False):
    """
    Ref.:
    -----------
    Schorpp, Kenji, et al.
    J Biomol Screen, 19.5 (2014): 715-726.
        
    Brief:
    -----------
    The endpoint HIS_FHs presents
    a compound prevents the binding of the protein His-tag moiety to nickel chelate
    There are 19 SMARTS in this endpoint
        
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_AlphaScreen_HIS_FHs(mols, detail=True)
    """
    HIS = _Filter('AlphaScreen_HIS_FHs',detail)
    HIS.get_pattl()
    res = list(map(HIS.scan,mols))
    return res


def Check_Biodegradable(mols, detail=False):
    """
    Ref:
    -----------
    Environment Canada.
    Existing Substances Program (CD-ROM), released April, 2004 (2003).
        
    Brief:
    -----------
    The endpoint Biodegradable presents
    a compound may be Biodegradable.
    There are 9 SMARTS in this enpoint
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Biodegradable(mols, detail=True)
    """
    Biodegradable = _Filter('Biodegradable',detail)
    Biodegradable.get_pattl()
    res = list(map(Biodegradable.scan,mols))
    return res
    

def Check_Chelating(mols, detail=False):
    """
    Ref.:
    -----------
    Agrawal, Arpita, et al.
    ChemMedChem, 5.2 (2010): 195-199.
    
    Brief:
    -----------
    The endpoint Chelating presents
    a compound may inhibit metalloproteins.
    Thers are 55 SMARTS in this endpoint
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Chelating(mols, detail=True)
    """
    Chelating = _Filter('Chelating',detail)
    Chelating.get_pattl()
    res = list(map(Chelating.scan,mols))
    return res


def Check_Developmental_Mitochondrial(mols, detail=False):
    """
    Ref.:
    -----------
    Mukesh, P.
    Structural Alerts for Developmental Toxicity and 
    Mitochondrial Toxicity Molecular Initiating Events (Lhasa Limited)
        
    Brief:
    -----------
    The endpoint Developmental_Mitochondrial present
    a compound may casue Developmental Toxicity and Mitochondrial Toxicity
    There are 12 SMARTS in this endpoint
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Developmental_Mitochondrial(mols, detail=True)
    """
    Developmental_Mitochondrial = _Filter('Developmental_Mitochondrial',detail)
    Developmental_Mitochondrial.get_pattl()
    res = list(map(Developmental_Mitochondrial.scan,mols))
    return res


def Check_Genotoxic_Carcinogenicity_Mutagenicity(mols, detail=False):
    """
    Ref.:
    -----------
    (1) Benigni, Romualdo, and Cecilia Bossa. 
        Mutat Res Rev Mutat Res, 659.3 (2008): 248-261.
    (2) Ashby, John, and Raymond W. Tennant.
        Mutat Res Rev Mutat Res, 204.1 (1988): 17-115.
    (3) Kazius, Jeroen, Ross McGuire, and Roberta Bursi.
        J Med Chem, 48.1 (2005): 312-320.
    (4) Bailey, Allan B., et al.
        Regul Toxicol Pharmacol, 42.2 (2005): 225-235.
    
    Brief:
    -----------
    The endpoint Genotoxic_Carcinogenicity_Mutagenicity present
    a compound may cause carcinogenicity or(and) mutagenicity through genotoxic mechanisms;
    There are 117 SMARTS in this endpoint
        
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Genotoxic_Carcinogenicity_Mutagenicity(mols, detail=True)
    """
    Genotoxic_Carcinogenicity_Mutagenicity = _Filter('Genotoxic_Carcinogenicity_Mutagenicity',detail)
    Genotoxic_Carcinogenicity_Mutagenicity.get_pattl()
    res = list(map(Genotoxic_Carcinogenicity_Mutagenicity.scan,mols))
    return res
    
    
def Check_Idiosyncratic(mols, detail=False):
    """
    Ref.:
    -----------
    Kalgutkar, Amit S., and John R. Soglia. 
    Expert Opin Drug Metab Toxicol, 1.1 (2005): 91-142
        
    Brief:
    -----------
    The endpoit Idiosyncratic presents
    a compound may has diosyncratic toxicity
    There are 35 SMARTS in this endpoint.
        
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Idiosyncratic(mols, detail=True)
    """
    Idiosyncratic = _Filter('Idiosyncratic',detail)
    Idiosyncratic.get_pattl()
    res = list(map(Idiosyncratic.scan,mols))
    return res


def Check_LD50_Oral(mol, detail=False):
    """
    Ref.:
    -----------
    Tinkov OV
    Biomed Khim, 65.2 (2019): 123-132.
       
    Brief:
    -----------
    The endpoint LD50_Oral presentshh
    a compound may cause acute toxicity during oral administration;
    There are 20 SMARTS in this endpoint
         
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_LD50_Oral(mols, detail=True)
    """
    LD50_Oral = _Filter('LD50_Oral',detail)
    LD50_Oral.get_pattl()
    res = list(map(LD50_Oral.scan,mols))
    return res
    

def Check_Luciferase_Inhibitory(mols, detail=False):
    """
    Ref.:
    -----------
    An unpublished parper
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Luciferase_Inhibitory(mols, detail=True)
    """
    Luciferase_Inhibitory = _Filter('Luciferase_Inhibitory',detail)
    Luciferase_Inhibitory.get_pattl()
    res = list(map(Luciferase_Inhibitory.scan,mols))
    return res
    

def Check_NonBiodegradable(mols, detail=False):
    """
    Ref:
    -----------
    Environment Canada.
    Existing Substances Program (CD-ROM), released April, 2004 (2003).
        
    Brief:
    -----------
    The endpoint Biodegradable presents
    a compound may be non-biodegradable.
    There are 19 SMARTS in this enpoint
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_NonBiodegradable(mols, detail=True)
    """
    NonBiodegradable = _Filter('NonBiodegradable',detail)
    NonBiodegradable.get_pattl()
    list(map(NonBiodegradable.scan,mols))
    return res
    

def Check_NonGenotoxic_Carcinogenicity(mols, detail=False):
    """
    Ref.:
    -----------
    (1) Benigni, Romualdo, and Cecilia Bossa. 
        Mutat Res Rev Mutat Res, 659.3 (2008): 248-261.
    (2) Benigni, Romualdo, Cecilia Bossa, and Olga Tcheremenskaia.
        Chem Rev, 113.5 (2013): 2940-2957.
    
    Brief:
    -----------
    The endpoint Genotoxic_Carcinogenicity_Mutagenicity present
    a compound may cause carcinogenicity or(and) mutagenicity through genotoxic mechanisms;
    There are 23 SMARTS in this endpoint
        
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_NonGenotoxic_Carcinogenicity(mols, detail=True)
    """
    NonGenotoxic_Carcinogenicity = _Filter('NonGenotoxic_Carcinogenicity',detail)
    NonGenotoxic_Carcinogenicity.get_pattl()
    res = list(map(NonGenotoxic_Carcinogenicity.scan,mols))
    return res


def Check_PAINS(mols, detail=False):
    """
    Ref.:
    -----------
    Baell, Jonathan B., and Georgina A. Holloway.
    J Med Chem, 53.7 (2010): 2719-2740.
    
    Brief:
    -----------
    PAINS i.e.,Pan Assay Interference Compounds, presents a type of
    compounds tend to be hitted in HTS. The article, J Med Chem, 53.7 (2010): 2719-2740,
    proposed 480 substructures to demonstrated PAINS.
           
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_PAINS(mols, detail=True)
    """
    PAINS = _Filter('Pains',detail)
    PAINS.get_pattl()
    res = list(map(PAINS.scan,mols))
    return res


def Check_Potential_Electrophilic(mols, detail=False):
    """
    Ref.:
    -----------
    Enoch, S. J., et al.
    Crit Rev Toxicol, 41.9 (2011): 783-802.
        
    Brief:
    -----------
    The point Potential_Electrophilic presents
    a compound would be more probably take part in electrophilic reaction,
    and the electrophilic reaction is strongly assosiated with protein binding.
    There are 119 SMARTS in this endpoint
        
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Potential_Electrophilic(mols, detail=True)
    """
    Potential_Electrophilic = _Filter('Potential_Electrophilic',detail)
    Potential_Electrophilic.get_pattl()
    res = list(map(Potential_Electrophilic.scan,mols))
    return res

def Check_Promiscuity(mols, detail=False):
    """
    Ref.:
    -----------
    Pearce, Bradley C., et al.
    J Chem Inf Model, 46.3 (2006): 1060-1068.
        
    Brief:
    -----------
    There are 119
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Promiscuity(mols, detail=True)
    """
    Promiscuity = _Filter('Promiscuity',detail)
    Promiscuity.get_pattl()
    res = list(map(Promiscuity.scan,mols))
    return res
    

def Check_Reactive_Unstable_Toxic(mols, detail=False):
    """
    Ref.:
    -----------
    ChemDiv
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Reactive_Unstable_Toxic(mols, detail=True)
    """
    Reactive_Unstable_Toxic = _Filter('Reactive_Unstable_Toxic',detail)
    Reactive_Unstable_Toxic.get_pattl()
    res = list(map(Reactive_Unstable_Toxic.scan,mols))
    return res


def Check_Skin_Sensitization(mols, detail=False):
    """
    Ref.:
    -----------
    (1) Payne, M. P., and P. T. Walsh
        J Chem Inf Comput Sci, 34.1 (1994): 154-161.
    (2) Enoch, S. J., J. C. Madden, and M. T. D. Cronin
        SAR QSAR Environ Res, 19.5-6 (2008): 555-578.
    (3) Barratt, M. D., et al.
        Toxicol In Vitro, 8.5 (1994): 1053-1060.
            
    Brief:
    -----------
    There are 155 SMARTS in this endpoint
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Reactive_Skin_Sensitization(mols, detail=True)
    """
    Skin_Sensitization = _Filter('Skin_Sensitization',detail)
    Skin_Sensitization.get_pattl()
    res = list(map(Skin_Sensitization.scan,mols))
    return res
    

def Check_DNA_Binding(mols, detail=False):
    """
    Ref.:
    -----------   

    Brief:
    -----------
    There are 78 SMARTS in this endpoint
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Reactive_DNA_Binding(mols, detail=True)
    """
    DNA_Binding = _Filter('DNA_Binding',detail)
    DNA_Binding.get_pattl()
    res = list(map(DNA_Binding.scan,mols))
    return res


def Check_SureChEMBL(mol, detail=False):
    """
    Ref.: https://www.surechembl.org/knowledgebase/169485   
    -----------
    
    Brief:
    -----------
    SMARTS descriptions used by SureChEMBL that would indicate whether 
    a compound would match one or more structural alerts and
    hence considered to have a MedChem unfriendly status.  
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_SureChEMBL(mols, detail=True)
    """
    SureChEMBL = _Filter('SureChEMBL',detail)
    SureChEMBL.get_pattl()
    res = list(map(SureChEMBL.scan,mols))
    return res
    
    
def Check_BMS(mols, detail=False):
    """
    Ref.:
    -----------
    Pearce, Bradley C., et al.
    J Chem Inf Model, 46.3 (2006): 1060-1068.
    
    Brief:
    -----------
    Pearce has proposed a Functional Group Compound Filters(FG Filters).
    The FG filters are consisted of two part, Exclusion FG filters and informational filters.
    Exclusion FG filters are those intended for compound removal from screening decks;
    Informational filters are useful for compound annotation.
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_BMS(mols, detail=True)
    """
    BMS = _Filter('BMS',detail)
    BMS.get_pattl()
    res = list(map(BMS.scan,mols))
    return res


def Check_NTD(mols, detail=False):
    """
    Ref.:
    -----------
    Brenk, Ruth, et al.
    ChemMedChem, 3(3), 435-444.
          
    Brief:
    -----------
    Brenk has proposed 105 unwanted groups in HTS
    
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_NTD(mols, detail=True)
    """        
    NTD = _Filter('NTD',detail)
    NTD.get_pattl()
    res = list(map(NTD.scan,mols))
    return res


def Check_Alarm_NMR(mols, detail=False):
    """
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Alarm_NMR(mols, detail=True)
    """
    Alarm_NMR = _Filter('Alarm_NMR',detail)
    Alarm_NMR.get_pattl()
    res = list(map(Alarm_NMR.scan,mols))
    return res
    

def Check_Frequent_Hitters(mols, detail=False):
    """
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Alarm_NMR(mols, detail=True)
    """
    Frequent_Hitters = _Filter('Frequent_Hitters',detail)
    Frequent_Hitters.get_pattl()
    res = list(map(Frequent_Hitters.scan,mols))
    return res
    

def Check_Aggregators(mols, detail=False):
    """
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Aggregators(mols, detail=True)
    """
    Aggregators = _Filter('Aggregators',detail)
    Aggregators.get_pattl()
    res = list(map(Aggregators.scan,mols))
    return res


def Check_Toxicophores(mols, stype='single', detail=False):
    """
    Parameters:
    -----------         
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
    -----------
    a lisf of namedtuple            
    namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
    and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    
    Usage:
    -----------
    smis = ['C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1']
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_Alarm_NMR(mols, detail=True)
    """
    Toxicophores = _Filter('Toxicophores',detail)
    Toxicophores.get_pattl()
    res = list(map(Toxicophores.scan,mols))
    return res
    

def VisualizeFragment(mol,highlightAtoms,figsize=[400,200]):
    from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D
    from IPython.display import SVG
    from rdkit.Chem import rdDepictor
    """
    This function is used for show which part of fragment matched the SMARTS
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
        the molecule to be visualized
    atoms: tuple
        the index of atoms to be highlighted
    
    Rrturn:
    -----------
    pic: IPython.core.display.SVG
        a SVG file
    
    Usage:
    ----------- 
    mol = Chem.MolFromSmiles('C1=CC=C2C(=O)CC(=O)C2=C1')    
    pic = VisualizeFragment(mol,(0, 1, 2, 6, 7, 8,10))
    """    
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(*figsize)
    drawer.DrawMolecule(mol,highlightAtoms=highlightAtoms)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    fig = SVG(svg)
    return fig



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
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    res = Check_PAINS(mols,detail=True)
    
    
    
    