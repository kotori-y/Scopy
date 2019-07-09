# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 21:59:42 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""


from rdkit.Chem import AllChem as Chem
from scopy.StructureAlert import SmartProcess


def Check_Acute_Aquatic_Toxicity(mol, detail=False):
    """
    ---
    Ref.:
        (1) Hermens, J. L.
        Environ Health Perspect, 87 (1990): 219-225.
        (2) Verhaar, Henk JM, Cees J. Van Leeuwen, and Joop LM Hermens.
        Chemosphere, 25.4 (1992): 471-491.
    
    ---
    Brief:
        The endpoint, Acute_Aquatic_Toxicity presents
        a compound may cause toxicity to liquid(water).
        There are 99 SMARTS in this endpoint.
   
    ---     
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple            
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Acute_Aquatic_Toxicity',detail=detail)    
    return res


def Check_AlphaScreen_FHs(mol, detail=False):
    """
    ---
    Ref.:
        Schorpp, Kenji, et al.
        J Biomol Screen, 19.5 (2014): 715-726.
        
    ---
    Brief:
        The endpoint, AlphaScreen_FHs presents
        a compound may be alphascreen frequent hitters
        There are 6 SMARTS in this endpoint
   
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='AlphaScreen_FHs',detail=detail) 
    return res


def Check_AlphaScreen_GST_FHs(mol, detail=False):
    """
    ---
    Ref.:
        Brenke, Jara K., et al.
        J Biomol Screen, 21.6 (2016): 596-607.
        
    ---
    Brief:
        The endpoint GST_FHs presents
        a compound may prevents GST/GSH interaction during HTS.
        There are 34 SMARTS in this endpoint
        
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='AlphaScreen_GST_FHs',detail=detail)  
    return res


def Check_AlphaScreen_HIS_FHs(mol, detail=False):
    """
    ---
    Ref.:
        Schorpp, Kenji, et al.
        J Biomol Screen, 19.5 (2014): 715-726.
        
    ---
    Brief:
        The endpoint HIS_FHs presents
        a compound prevents the binding of the protein His-tag moiety to nickel chelate
        There are 19 SMARTS in this endpoint
        
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='AlphaScreen_HIS_FHs') 
    return res


def Check_Biodegradable(mol, detail=False):
    """
    ---
    Ref:
        Environment Canada.
        Existing Substances Program (CD-ROM), released April, 2004 (2003).
        
    ---
    Brief:
        The endpoint Biodegradable presents
        a compound may be Biodegradable.
        There are 9 SMARTS in this enpoint
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Biodegradable',detail=detail)    
    return res
    

def Check_Chelating(mol, detail=False):
    """
    ---
    Ref.:
        Agrawal, Arpita, et al.
        ChemMedChem, 5.2 (2010): 195-199.
    
    ---
    Brief:
        The endpoint Chelating presents
        a compound may inhibits metalloproteins.
        Thers are 55 SMARTS in this endpoint
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Chelating',detail=detail)
    return res


def Check_Developmental_Mitochondrial(mol, detail=False):
    """
    Ref.:
        Mukesh, P.
        Structural Alerts for Developmental Toxicity and 
        Mitochondrial Toxicity Molecular Initiating Events (Lhasa Limited)
        
    ---
    Brief:
        The endpoint Developmental_Mitochondrial present
        a compound may casue Developmental Toxicity and Mitochondrial Toxicity
        There are 12 SMARTS in this endpoint
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Developmental_Mitochondrial')
    return res


def Check_Genotoxic_Carcinogenicity_Mutagenicity(mol, detail=False):
    """
    ---
    Ref.:
        (1) Benigni, Romualdo, and Cecilia Bossa. 
            Mutat Res Rev Mutat Res, 659.3 (2008): 248-261.
        (2) Ashby, John, and Raymond W. Tennant.
            Mutat Res Rev Mutat Res, 204.1 (1988): 17-115.
        (3) Kazius, Jeroen, Ross McGuire, and Roberta Bursi.
            J Med Chem, 48.1 (2005): 312-320.
        (4) Bailey, Allan B., et al.
            Regul Toxicol Pharmacol, 42.2 (2005): 225-235.
    
    ---
    Brief:
        The endpoint Genotoxic_Carcinogenicity_Mutagenicity present
        a compound may cause carcinogenicity or(and) mutagenicity through genotoxic mechanisms;
        There are 117 SMARTS in this endpoint
        
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Genotoxic_Carcinogenicity_Mutagenicity',detail=detail)   
    return res
    
    
def Check_Idiosyncratic(mol, detail=False):
    """
    ---
    Ref.:
        Kalgutkar, Amit S., and John R. Soglia. 
        Expert Opin Drug Metab Toxicol, 1.1 (2005): 91-142
        
    ---
    Brief:
        The endpoit Idiosyncratic presents
        a compound may has diosyncratic toxicity
        There are 35 SMARTS in this endpoint.
        
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Idiosyncratic',detail=detail)
    return res
    

def Check_LD50_Oral(mol, detail=False):
    """
    ---
    Ref.:
        Tinkov OV
        Biomed Khim, 65.2 (2019): 123-132.
    
    ---    
    Brief:
        The endpoint LD50_Oral presents
         a compound may cause acute toxicity during oral administration;
         There are 20 SMARTS in this endpoint
         
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='LD50_Oral',detail=detail)
    return res


def Check_Luciferase_Inhibitory(mol, detail=False):
    """
    ---
    Ref.:
        An unpublished parper
    
    ---
    Brief:
        
        
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Luciferase_Inhibitory')
    return res
    

def Check_NonBiodegradable(mol, detail=False):
    """
    ---
    Ref:
        Environment Canada.
        Existing Substances Program (CD-ROM), released April, 2004 (2003).
        
    ---
    Brief:
        The endpoint Biodegradable presents
        a compound may be non-biodegradable.
        There are 19 SMARTS in this enpoint
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='NonBiodegradable',detail=detail)   
    return res
    
    
def Check_NonGenotoxic_Carcinogenicity(mol, detail=False):
    """
    ---
    Ref.:
        (1) Benigni, Romualdo, and Cecilia Bossa. 
            Mutat Res Rev Mutat Res, 659.3 (2008): 248-261.
        (2) Benigni, Romualdo, Cecilia Bossa, and Olga Tcheremenskaia.
            Chem Rev, 113.5 (2013): 2940-2957.
    
    ---
    Brief:
        The endpoint Genotoxic_Carcinogenicity_Mutagenicity present
        a compound may cause carcinogenicity or(and) mutagenicity through genotoxic mechanisms;
        There are 23 SMARTS in this endpoint
        
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='NonGenotoxic_Carcinogenicity',detail=detail)
    return res


def Check_Pains(mol, detail=False):
    """
    ---
    Ref.:
        Baell, Jonathan B., and Georgina A. Holloway.
        J Med Chem, 53.7 (2010): 2719-2740.
    
    ---
    Brief: 
        PAINS i.e.,Pan Assay Interference Compounds, presents a type of
        compounds tend to be hitted in HTS. The article, J Med Chem, 53.7 (2010): 2719-2740,
        proposed 480 substructures to demonstrated PAINS.
           
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False   
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Pains',detail=detail)
    return res


def Check_Potential_Electrophilic(mol, detail=False):
    """
    ---
    Ref.:
        Enoch, S. J., et al.
        Crit Rev Toxicol, 41.9 (2011): 783-802.
        
    ---
    Brief:
        The point Potential_Electrophilic presents
        a compound would be more probably take part in electrophilic reaction,
        and the electrophilic reaction is strongly assosiated with protein binding.
        There are 119 SMARTS in this endpoint
        
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False 
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Potential_Electrophilic',detail=detail)
    return res


def Check_Promiscuity(mol, detail=False):
    """
    ---
    Ref.:
        Pearce, Bradley C., et al.
        J Chem Inf Model, 46.3 (2006): 1060-1068.
        
    ---
    Brief:
        There are 119
    
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Promiscuity',detail=detail)
    
    return df


def Check_Reactive_Unstable_Toxic(mol, detail=False):
    """
    Ref.:
        ChemDiv
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Reactive_Unstable_Toxic',detail=detail)
    
    return df


def Check_Skin_Sensitization(mol, detail=False):
    """
    ---
    Ref.:
        (1) Payne, M. P., and P. T. Walsh
            J Chem Inf Comput Sci, 34.1 (1994): 154-161.
        (2) Enoch, S. J., J. C. Madden, and M. T. D. Cronin
            SAR QSAR Environ Res, 19.5-6 (2008): 555-578.
        (3) Barratt, M. D., et al.
            Toxicol In Vitro, 8.5 (1994): 1053-1060.
            
    ---
    Brief:
        There are 155 SMARTS in this endpoint
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Skin_Sensitization',detail=detail)    
    return res
    

def Check_DNA_Binding(mol,detail=False):
    """
    ---
    Ref.:
        
    ---
    Brief:
        There are 78 SMARTS in this endpoint
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False  
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='DNA_Binding',detail=detail)
    return res
    

def Check_SureChEMBL(mol,detail=False):
    """
    ---
    Ref.: https://www.surechembl.org/knowledgebase/169485   
    
    ---
    Brief: SMARTS descriptions used by SureChEMBL that would indicate whether 
    a compound would match one or more structural alerts and
    hence considered to have a MedChem unfriendly status.  
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='SureChEMBL',detail=detail)
    return res
    
    
def Check_BMS(mol,detail=False):
    """
    ---
    Ref.: Pearce, Bradley C., et al.
          J Chem Inf Model, 46.3 (2006): 1060-1068.
    
    ---
    Brief: Pearce has proposed a Functional Group Compound Filters(FG Filters).
    The FG filters are consisted of two part, Exclusion FG filters and informational filters.
    Exclusion FG filters are those intended for compound removal from screening decks;
    Informational filters are useful for compound annotation.
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='BMS',detail=detail)
    return res


def Check_NTD(mol,detail=False):
    """
    ---
    Ref.: Brenk, Ruth, et al.
          ChemMedChem, 3(3), 435-444.
          
    ---
    Brief: Brenk has proposed 105 unwanted groups in HTS
    
    ---
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
    
    Return:
        a namedtuple
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """        
    res = SmartProcess._CheckWithSmarts(mol,endpoint='NTD',detail=detail)
    return res


if __name__ == '__main__':
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
    
    smi = smis[-1]
    mol = Chem.MolFromSmiles(smi)
    res = Check_AlphaScreen_GST_FHs(mol,detail=True)
    #res = Check_Acute_Aquatic_Toxicity(mol)
    print(res)   