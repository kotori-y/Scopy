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
    #################################################################
    This function is used for detecting molecular whether or not
    has the toxicophores of Acute Aquatic Toxicity through SMARTS.
            
    -Endpoint Info:        
        
        Endpoint: Acute Aquatic Toxicity
        
        Brief Express: The toxicity to water
                     
    -Ref.:
        
        (1) Hermens, J. L.
        Environmental health perspectives 87 (1990): 219-225.
        (2) Verhaar, Henk JM, Cees J. Van Leeuwen, and Joop LM Hermens.
        Chemosphere 25.4 (1992): 471-491.
        
    -Usage:
        
        res = Check_Acute_Aquatic_Toxicity(mol)
        
        Input: mol is a molecule object.
        
        Output: res is a customize object--CheckRes            
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Acute_Aquatic_Toxicity',detail=detail)    
    return res


def Check_AlphaScreen_FHs(mol, detail=False):
    """
    #################################################################
    This Function is used for detecting molecular whether or not
    has the substructure of alphascreen frequent hitters, through SMARTS,
        
    We have obtained 6 substructures in this endpoint
    
    -Ref.:
        Schorpp, Kenji, et al.
        Journal of biomolecular screening 19.5 (2014): 715-726.
        
    -Usage:
        
        df = Check_AlphaScreen_FHs(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='AlphaScreen_FHs',detail=detail)
    
    return res


def Check_AlphaScreen_GST_FHs(mol, detail=False):
    """
    NOTE: The Smarts in this endpoint may should be revised
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='AlphaScreen_GST_FHs',detail=detail)
    
    return res


def Check_AlphaScreen_HIS_FHs(mol, detail=False):
    """
    #################################################################
    This Function is used for detecting molecular whether or not
    has the substructure of alphascreen frequent hitters in HIS, through SMARTS,
        
    We have obtained 19 substructures in this endpoint
    
    -Ref.:
        Schorpp, Kenji, et al.
        Journal of biomolecular screening 19.5 (2014): 715-726.
        
    -Usage:
        
        df = Check_AlphaScreen_HIS_FHs(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='AlphaScreen_HIS_FHs')
    
    return res


def Check_Biodegradable(mol, detail=False):
    """
    #################################################################
    This function is used for detecting molecular whether or not
    be able to be biodegraded
    
    We have obtained 9 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Biodegradable(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Biodegradable',detail=detail)
    
    return df
    

def Check_Chelating(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 55 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Chelating(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Chelating',detail=detail)
    
    return res

def Check_Developmental_Mitochondrial(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 12 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Developmental_Mitochondrial(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Developmental_Mitochondrial')

    return df


def Check_Genotoxic_Carcinogenicity_Mutagenicity(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 117 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Genotoxic_Carcinogenicity_Mutagenicity(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Genotoxic_Carcinogenicity_Mutagenicity',detail=detail)
    
    return df
    
    
def Check_Idiosyncratic(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 35 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Idiosyncratic(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Idiosyncratic',detail=detail)
    
    return res
    

def Check_LD50_Oral(mol, detail=False):
    """
    #################################################################
    This function is used for...    
    
    We have obtained 20 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_LD50_oral(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='LD50_Oral',detail=detail)
    
    return res


def Check_Luciferase_Inhibitory(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 3 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Luciferase_Inhibitory(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Luciferase_Inhibitory')
    
    return df
    

def Check_NonBiodegradable(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 19 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_NonBiodegradable(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='NonBiodegradable',detail=detail)
    
    return res
    
    
def Check_NonGenotoxic_Carcinogenicity(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 23 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_NonGenotoxic_Carcinogenicity(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='NonGenotoxic_Carcinogenicity',detail=detail)
    
    return df


def Check_Pains(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 480 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        res = Check_Pains(mol)
        
        Input: mol is a molecule object.
        
        Output: res is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Pains',detail=detail)
    
    return res


def Check_Potential_Electrophilic(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 119 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Potential_Electrophilic(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Potential_Electrophilic',detail=detail)
    
    return res

def Check_Promiscuity(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained  substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Potential_Electrophilic(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Promiscuity',detail=detail)
    
    return df


def Check_Reactive_Unstable_Toxic(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 335 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Reactive_Unstable_Toxic(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Reactive_Unstable_Toxic',detail=detail)
    
    return df


def Check_Skin_Sensitization(mol, detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 155 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_Reactive_Unstable_Toxic(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    df = SmartProcess._CheckWithSmarts(mol,endpoint='Skin_Sensitization',detail=detail)    
    
    return df
    

def Check_DNA_Binding_ToxTree(mol,detail=False):
    """
    #################################################################
    This function is used for...
    
    We have obtained 78 substructures in this endpoint
    
    -Ref.:
        
    -Usage:
        
        df = Check_DNA_Binding_ToxTree(mol)
        
        Input: mol is a molecule object.
        
        Output: df is a pandas.core.frame.DataFrame object.
    #################################################################
    """
    res = SmartProcess._CheckWithSmarts(mol,endpoint='Acute_Aquatic_Toxicity',detail=detail)
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
    
    Return: a namedtuple
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
    
    Return: a namedtuple
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
        >>> mol: rdkit.Chem.rdchem.Mol.
        >>> detail: bool, optional(default=True), 
            When set to True, function will return more information(MatchedAtoms,MatchedNames)
            else, only return Disposed and Endpoint
    
    Return: a namedtuple
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
    res = Check_SureChEMBL(mol,detail=True)
    #res = Check_Acute_Aquatic_Toxicity(mol)
    print(res)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    