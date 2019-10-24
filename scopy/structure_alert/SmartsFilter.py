# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:48:23 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

try:
	from . import FilterWithSmarts
except:
	import FilterWithSmarts
from rdkit.Chem import AllChem as Chem


class Filter(object):
    """
    Here, we comdat the whole function that check endpoint retrieved  from module FilterWithSmarts
    
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    """
    def __init__(self, mols):
        self.mols = mols
        
    def Check_Acute_Aquatic_Toxicity(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Acute_Aquatic_Toxicity = FilterWithSmarts.Check_Acute_Aquatic_Toxicity(self.mols, detail)
        
    def Check_Check_AlphaScreen_FHs(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.AlphaScreen_FHs = FilterWithSmarts.Check_AlphaScreen_FHs(self.mols, detail)
        
    def Check_AlphaScreen_GST_FHs(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.AlphaScreen_GST_FHs = FilterWithSmarts.Check_AlphaScreen_GST_FHs(self.mols, detail)
    
    def Check_AlphaScreen_HIS_FHs(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.AlphaScreen_HIS_FHs = FilterWithSmarts.Check_AlphaScreen_HIS_FHs(self.mols, detail)
        
    def Check_Biodegradable(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Biodegradable = FilterWithSmarts.Check_Biodegradable(self.mols, detail)
        
    def Check_Chelating(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Chelating = FilterWithSmarts.Check_Chelating(self.mols, detail)
        
    def Check_Developmental_Mitochondrial(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Developmental_Mitochondrial = FilterWithSmarts.Check_Developmental_Mitochondrial(self.mols, detail)
        
    def Check_Genotoxic_Carcinogenicity_Mutagenicity(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Genotoxic_Carcinogenicity_Mutagenicity = FilterWithSmarts.Check_Genotoxic_Carcinogenicity_Mutagenicity(self.mols, detail)
        
    def Check_Idiosyncratic(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Check_Idiosyncratic = FilterWithSmarts.Check_Idiosyncratic(self.mols, detail)
        
    def Check_LD50_Oral(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Check_LD50_Oral = FilterWithSmarts.Check_LD50_Oral(self.mols, detail)
        
    def Check_Luciferase_Inhibitory(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Luciferase_Inhibitory = FilterWithSmarts.Check_Luciferase_Inhibitory(self.mols, detail)
        
    def Check_NonBiodegradable(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.NonBiodegradable = FilterWithSmarts.Check_NonBiodegradable(self.mols, detail)
        
    def Check_NonGenotoxic_Carcinogenicity(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.NonGenotoxic_Carcinogenicity = FilterWithSmarts.Check_NonGenotoxic_Carcinogenicity(self.mols, detail)
        
    def Check_PAINS(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.PAINS = FilterWithSmarts.Check_PAINS(self.mols, detail)
        
    def Check_Potential_Electrophilic(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Potential_Electrophilic = FilterWithSmarts.Check_PAINS(self.mols, detail)
        
    def Check_Promiscuity(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Promiscuity = FilterWithSmarts.Check_Promiscuity(self.mols, detail)
    
    def Check_Reactive_Unstable_Toxic(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Reactive_Unstable_Toxic = FilterWithSmarts.Check_Reactive_Unstable_Toxic(self.mols, detail)
    
    def Check_Skin_Sensitization(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Skin_Sensitization = FilterWithSmarts.Check_Skin_Sensitization(self.mols, detail)
        
    def Check_DNA_Binding(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.DNA_Binding = FilterWithSmarts.Check_DNA_Binding(self.mols, detail)
        
    def Check_SureChEMBL(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.SureChEMBL = FilterWithSmarts.Check_SureChEMBL(self, detail)
       
    def Check_BMS(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.BMS = FilterWithSmarts.Check_BMS(self.mols, detail)
        
    def Check_NTD(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.NTD = FilterWithSmarts.Check_NTD(self.mols, detail)
    
    def Check_Alarm_NMR(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Alarm_NMR = FilterWithSmarts.Check_Alarm_NMR(self.mols, detail)
    
    def Check_Frequent_Hitters(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Frequent_Hitters = FilterWithSmarts.Check_Frequent_Hitters(self.mols, detail)

    def Check_Aggregators(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Aggregators = FilterWithSmarts.Check_Aggregators(self.mols, detail)
     
    def Check_Toxicophores(self, detail=False):
        """
        Parameters:
        -----------
        detail: bool, optional(default=True), 
        When set to True, function will return more information(MatchedAtoms,MatchedNames)
        else, only return Disposed and Endpoint
        """
        self.Toxicophores = FilterWithSmarts.Check_Toxicophores(self.mols, detail)
        
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
    F = Filter(mols)
    F.Check_Aggregators()
    F.Check_Alarm_NMR()
    print(F.Aggregators)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    