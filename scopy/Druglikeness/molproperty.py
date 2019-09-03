# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 21:59:42 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

__doc__ = """
    This moudle is used to calculated properise that contained in our collectded rules
    ---
    up to now(2019.07.22), we have achived followed properties:
        Molcular Weight >>> MW
        Number of bonds >>> nBond
        Number of atoms >>> nAtom
        Number of heteroatoms >>> nHet
        Number of heavy atom >>> nHev
        Number of rotable bonds >>> nRot
        Number of rigid bonds >>> nRig
        Number of SSSR >>> nRing
        logP >>> logP
        logD >>> logD
        logSw >>> logSw
        Acid or Base >>> ab
        pKa >>> pKa
        QED >>> qed
        Molecular refraction >>> MR
        Number of hydrogen bond donors >>> nHD
        Number of hydrogen bond acceptors >>> nHA
        Number of hydrogen bond donors& acceptors >>> nHB
        Aromatic proportion >>> AP
        sp3 hybridized carbons/total carbon count >>> Fsp3
        TPSA >>> TPSA
        Number of atoms involved in the biggest system ring >>> MaxRing
        Number of Sterocenterss >>> nStero
        HetCarbonRatio >>> HetRatio
        synthetic accessibility score >>> SAscore
        natural product-likeness score >>> NPscore
        Number of single bonds >>> nSingle
        Number of double bobds >>> nDoudle
        Number of triple bonds >>> nTriple
        Volume of mol >>> Vol
        Density >>> Dense
        MolFCharge >>> fChar
        Number of Carbon atoms >>> nC
        Number of Boron atoms >>> nB
        Number of Chlorin atoms >>> nCl
        Number of Bromine atoms >>> nBr
        Number of Iodine atoms >>> nI
        Number of Phosphor atoms >>> P
        Number of Sulfur atoms >>> nS
        Number of Oxygen atoms >>> nO
        Number of Nitrogen atoms >>> nN   
    ---
    Followed should be achieved in the future:
        logD
        difference between clogP and alogP
        Formal total charge of the compound
        Number of Charged Groups
    """
import sys,os
sys.path.append('..')
from fingerprint.fingerprints import CalculateGhoseCrippen
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors, Lipinski, QED
from itertools import combinations
from collections import namedtuple
from rdkit import RDConfig
ContriDir = RDConfig.RDContribDir
import csv
sys.path.append(ContriDir)
import ScoConfig
from SA_Score import sascorer
from NP_Score import npscorer
from IFG.ifg import identify_functional_groups



def CalculateMolWeight(mol):    
    """
    #################################################################
    Calculation of molecular weight
        
    ---->MW  
    
    Usage:
        
        result=CalculateMolWeight(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
        
    #################################################################
    """
    mol = Chem.AddHs(mol)
    MW = round(sum([atom.GetMass() for atom in mol.GetAtoms()]),2)
    return MW


def CalculateNumBonds(mol):
    """
    #################################################################
    Calculation the number of bonds where between heavy atoms
       
    --->nBond
    
    Usage:
        
        result = CalculateNumBonds(mol)
        
        Input: mol is a molecular object.
        
        Output: result is a numeric value.
    #################################################################
    """
    nBond = mol.GetNumBonds()    
    return nBond


def CalculateNumAtoms(mol):
    """
    #################################################################
    Calculation of the number of atoms in molecular
    
    ---->nAtom
    
    Usage:
        
        result = CalculateNumAtom(mol)
        
        Input: mol is a molecular object.
        
        Output: result is a numeric value.
    #################################################################
    """  
    mol = Chem.AddHs(mol)
    return mol.GetNumAtoms()
   

def CalculateHeteroNumber(mol):
    """
    #################################################################
    Calculation of Hetero counts in a molecule
    
    ---->nHet
    
    Usage:
        
        result=CalculateHeteroNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    i = len(
            [atom for atom in mol.GetAtoms()\
             if atom.GetAtomicNum()==6 or atom.GetAtomicNum()==1]
            )
    return mol.GetNumAtoms()-i


def CalculateNumRotatableBonds(mol):
    """
    #################################################################
    Calculation of the number of rotatableBonds
    
    NOTE: In some situaion Amide C-N bonds are not considered 
          because of their high rotational energy barrier
    
    ---->nRot
    
    Usage:
        
        result = CalculateNumRotatableBonds(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric values
    #################################################################
    """
    patt = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    nROT = len(mol.GetSubstructMatches(patt))    
    return nROT


def CalculateNumRigidBonds(mol):
    """
    #################################################################
    Number of non-flexible bonds, in opposite to rotatable bonds
    
    NOTE: This function need to be revised in the future(2019/05/09) 
    
    ---->nRig
    
     Usage:
        
        result = CalculateNumRigidBonds(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric value
    #################################################################
    """
    nBOND = CalculateNumBonds(mol)    
    flex = 0
    for bond in mol.GetBonds():
        bondtype = bond.GetBondType()
        if bondtype == Chem.rdchem.BondType.SINGLE and not bond.IsInRing():
            flex+=1               
    return nBOND-flex   


def CalculateNumRing(mol):
    """
    #################################################################
    Calculation of the number of ring
    
    ---->nRing
    
    Usage:
        
        result = CalculateRingNum(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric value
    #################################################################
    """
    nRing = Chem.GetSSSR(mol)    
    return nRing


def CalculateNumHeavyAtom(mol):
    """
    #################################################################
    Calculation of Heavy atom counts in a molecule
    
    ---->nHev
    
    Usage:
        
        result=CalculateHeavyAtomNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return mol.GetNumHeavyAtoms()


def CalculateLogD(mol):
    """
    #################################################################
    Calculation of molecular logD
    
    ---->LogD
    
    Usage:
        
        result = CalculateLogD(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric value
    #################################################################
    """
    intercept = 0.5748907159915493
    
    fps = CalculateGhoseCrippen([mol]).flatten()
    with open(ScoConfig.CrippenDir + '\\Crippen.txt') as f_obj:
        lines = csv.reader(f_obj,delimiter='\t')
        next(lines)
        contri = [x[-1] for x in lines]
        contri = [float(x) for x in contri]
    f_obj.close()
    logD = sum([a*b for a,b in zip(fps,contri)]) + intercept
    return logD


def CalculateLogP(mol):    
    """
    #################################################################
    Calculation of molecular LogP
    
    ---->logP
    
    Usage:
        
        result = CalculateLogP(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric value
    #################################################################
    """   
    return round(Descriptors.MolLogP(mol),2)  


def CheckAcid(mol):
    """
    """
    acid_fragment = ['[!H0;F,Cl,Br,I,N+,$([OH]-*=[!#6]),+]',
                     '[CX3](=O)[OX2H1]',
                     '[CX3](=O)[OX1H0-,OX2H1]',
                     '[$([OH]-*=[!#6])]',
                     '[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]',
                     '[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]',
                     '[CX3](=[OX1])[F,Cl,Br,I]'
                     ]
    for sma in acid_fragment:
        patt = Chem.MolFromSmarts(sma)
        if mol.HasSubstructMatch(patt):
            return 'acid'
    else:
        return 'base'        
        
def CalculatepKa(mol):
    from math import log10
    """
    Eq.:
        |pH-pKa| = log10(10^(logP-logD)-1)
        pKa = pH - log10(10^(logP-logD)-1) for acid
        pKa = log10(10^(logP-logD)-1) - pH for base  
    """
    logP = CalculateLogP(mol)
    logD = CalculateLogD(mol)
    status = CheckAcid(mol)
    try:
        if status == 'acid':
            pKa = 7.4 - log10(10**(logP-logD)-1)
        else:
            pKa = log10(10**(logP-logD)-1) - 7.4
        return pKa
    except:
        return 'N/A'
    
    
def CalculateMolMR(mol):
    """
    #################################################################
    Cacluation of molecular refraction value based on Crippen method
    
    ---->MR
    
    Usage:
        
        result = CalculateMolMR(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return round(Descriptors.MolMR(mol),2) 


def CalculateNumHDonors(mol):    
    """
    #################################################################
    Caculation of the number of Hydrogen Bond Donors
    
    ---->nHD
    
    Usage:
        
        result = CaculateNumHDonors(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric value
    #################################################################
    """
    nHD = Lipinski.NumHDonors(mol)    
    return nHD


def CalculateNumHAcceptors(mol):    
    """
    #################################################################
    Caculation of the number of Hydrogen Bond Acceptors
    
    ---->nHA
    
    Usage:
        
        result = CalculateNumHAcceptors(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric values
    #################################################################    
    """
    nHA = Lipinski.NumHAcceptors(mol)    
    return nHA


def CalculateNumHyBond(mol):
    """
    #################################################################
    Sum of Hydrogen Bond Donnors and Acceptors
    
    ---->nHB
    
     Usage:
        
        result = CaculateNumHyBond(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric value
    #################################################################
    """
    nHD = CalculateNumHDonors(mol)
    nHA = CalculateNumHAcceptors(mol)
    nHB = nHD+nHA
    return nHB


def CalculateAromaticProportion(mol):
    """
    #################################################################
    The proportion of heavy atoms in the molecule that are in an aromatic ring
    
    ---->AP
    
    Usage:
        
        result = CalculateAromaticProportion(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric value
    #################################################################
    """
    aroma = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a]')))
    total = CalculateNumHeavyAtom(mol)
    return round(aroma/total,2)    
    

def CalculateLogSw(mol):
    """
    #################################################################
    The logSw represents the logarithm of compounds water solubility computed by the ESOL method
    
    Equation: Log(Sw) = 0.16-0.638*clogP-0.0062*MWT+0.066*RB-0.74*AP
              >MWT: Molecular Weight
              >RB: Rotatable bonds
              >AP: Aromatic proportion
              
    Ref.: 
        Delaney, John S. 
        Journal of chemical information and computer sciences 44.3 (2004): 1000-1005.
    #################################################################
    """
    #Calculate each property
    MWT = CalculateMolWeight(mol)
    RB = CalculateNumRotatableBonds(mol)
    AP = CalculateAromaticProportion(mol)
    logP = CalculateLogP(mol)  
    logSw = 0.16-0.638*logP-0.0062*MWT+0.066*RB-0.74*AP   
    return round(logSw,2)


def CalculateFsp3(mol):
    """
    #################################################################
    Fsp3 (carbon bond saturation) is defined as the number of sp3 hybridized carbons / total carbon count.
    
    ---->FSP3
    
    Usage:
        
        result = CalculateFsp3(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return round(Lipinski.FractionCSP3(mol),2)
    

def CalculateTPSA(mol):
    """
    #################################################################
    Calculation of TPSA
    
    ---->tPSA
    
    Usage:
        
        result = CalculateTpsa(mol)
        
        Input: mol is a molecular object
        
        Output: result is a nmeric value
    #################################################################
    """
    return round(Descriptors.TPSA(mol),2)
    
            
def CalculateQED(mol,wtype = "mean"):
    """
    #################################################################
    Calculation QED descriptor under different weights
    
    A descriptor a measure of drug-likeness based on the concept of desirability
    -Ref.: Bickerton, G. Richard, et al.
           Nat Chem, 4.2 (2012): 90.
            
    Quantitative Estimate of Drug-likeness 
    
    ---->qed
    
    Usage:
        
        result = CalculateQED(mol,wtype='mean')
        
        Input: mol is a molecular object
        
        Output: result is a numeric values
    #################################################################
    """
    if wtype == "mean":
        qed = QED.weights_mean(mol)
    elif wtype == "max":
        qed = QED.weights_max(mol)
    elif wtype == "none":
        qed = QED.weights_none(mol)
    else:
        #msg = "invalid wtype has been input"
        qed = None
    
    return round(qed,2) 
    

def CalculateMaxSizeSystemRing(mol):
    """
    #################################################################
    Number of atoms involved in the biggest system ring
    
    ----> maxring
    
    Usage:
        
        result = CalculateMaxSizeSystemRing(mol)
        
        Input: mol is a molecular object
        
        Output: result is a numeric values
    #################################################################
    """
    #1.Obtaining which atoms consist of rings
    maxring = 0
    ri = mol.GetRingInfo()
    atoms = list(ri.AtomRings())    
    length = len(atoms)    
    if length == 0:
        pass
    else:
        rw = Chem.RWMol(mol)        
        #2.Judge which atoms are replacement
        atoms = [set(x) for x in atoms]            
        for pair in combinations(range(length),2):
            replace = list(atoms[pair[0]]&atoms[pair[1]])
            if len(replace) >= 2:
                for repl in list(combinations(replace,2)):
                    rw.RemoveBond(*repl)
            else:
                pass    
        m = Chem.MolFromSmiles(Chem.MolToSmiles(rw))
        ri = m.GetRingInfo()
        bonds = ri.BondRings()
        for item in bonds:
            if len(item) > maxring:
                maxring = len(item)   
    return maxring
    

def CalculateNumberStereocenters(mol):
    """
    #################################################################
    ----> nSTERO
    #################################################################
    """
    return Chem.CalcNumAtomStereoCenters(mol)    


def CalculateSolubilityForecastIndex(mol):
    """
    #################################################################
    Water Solubility estimated by Hill et al method
    
    ----> SFI
    
    SFI = clogD(pH=7.4) + #Ar
    #################################################################
    """
    logD = CalculateLogD(mol)
    pass


def _CalculateElementNumber(mol,AtomicNumber=6):
    """
    #################################################################
    **Internal used only**
    
    Calculation of element counts with atomic number equal to n in a molecule
    #################################################################
    """
    return len(
            [atom for atom in mol.GetAtoms()\
                if atom.GetAtomicNum() == AtomicNumber]
            )


def CalculateCarbonNumber(mol):
    """
    #################################################################
    Calculation of Carbon number in a molecule
    
    ---->ncarb
    
    Usage:
        
        result = CalculateCarbonNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """

    return _CalculateElementNumber(mol,AtomicNumber=6)


def CalculateBoronNumber(mol):

    """
    #################################################################
    Calculation of Boron counts in a molecule
    
    ---->ncof
    
    Usage:
        
        result=CalculateBoronNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
            
    return _CalculateElementNumber(mol,AtomicNumber=5)


def CalculateFluorinNumber(mol):

    """
    #################################################################
    Calculation of Fluorin counts in a molecule
    
    ---->ncof
    
    Usage:
        
        result=CalculateBoronNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
            
    return _CalculateElementNumber(mol,AtomicNumber=10)


def CalculateChlorinNumber(mol):

    """
    #################################################################
    Calculation of Chlorin counts in a molecule
    
    ---->ncocl
    
    Usage:
        
        result=CalculateChlorinNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """

    return _CalculateElementNumber(mol,AtomicNumber=17)


def CalculateBromineNumber(mol):

    """
    #################################################################
    Calculation of Bromine counts in a molecule
    
    ---->ncobr
    
    Usage:
        
        result=CalculateBromineNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """

    return _CalculateElementNumber(mol,AtomicNumber=35)


def CalculateIodineNumber(mol):
    """
    #################################################################
    Calculation of Iodine counts in a molecule
    
    ---->ncoi
    
    Usage:
        
        result=CalculateIodineNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """

    return _CalculateElementNumber(mol,AtomicNumber=53)


def CalculatePhosphorNumber(mol):
    """
    #################################################################
    Calcualtion of Phosphor number in a molecule
    
    ---->nphos
    
    Usage:
        
        result=CalculatePhosphorNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementNumber(mol,AtomicNumber=15)


def CalculateSulfurNumber(mol):
    """
    #################################################################
    Calculation of Sulfur counts in a molecule
    
    ---->nsulph
    
    Usage:
        
        result=CalculateSulfurNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementNumber(mol,AtomicNumber=16)


def CalculateOxygenNumber(mol):
    """
    #################################################################
    Calculation of Oxygen counts in a molecule
    
    ---->noxy
    
    Usage:
        
        result=CalculateOxygenNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################

    """
    return _CalculateElementNumber(mol,AtomicNumber=8)
        

def CalculateNitrogenNumber(mol):
    """
    #################################################################
    Calculation of Nitrogen counts in a molecule
    
    ---->nnitro
    
    Usage:
        
        result=CalculateNitrogenNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    
    return _CalculateElementNumber(mol,AtomicNumber=7)


def CalculateNumberChargedGroups(mol):
    """
    #################################################################
    Number of Charged Groups
    
    ---->ncharged
    
    Usage:
        
        result = CalculateNumberChargedGroups(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    pass
    
    return None


def CalculateTotalCharge(mol):
    """
    #################################################################
    Formal total charge of the compound.
    
    ---->totalchar
    
    Usage:
        
        result = CalculateTotalCharge(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    pass    
    return None


def CalculateHetCarbonRatio(mol):
    """
    #################################################################
    The ratio between the number of non carbon atoms and the number of carbon atoms.
    
    ---->HetCar
    
    Usage:
        
        result = CalculateHetCarbonRatio(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    Total = CalculateNumHeavyAtom(mol)
    nCarb = CalculateCarbonNumber(mol)
    het = Total-nCarb  
    return round(het/nCarb,2)    
    

def CalculateSAscore(mol):
    """
    ---
    Ref.: Ertl, Peter, and Ansgar Schuffenhauer.
          J Cheminform, 1(1), 8.      
    Based: https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score   
    
    ---
    Brief: A function to estimate ease of synthesis (synthetic accessibility) of drug-like molecules  
    
    ---
    Parameters:
        >>> mol: dkit.Chem.rdchem.Mol;
    
    Return:
        float, meant SA score
    """
    return round(sascorer.calculateScore(mol),2)


def CalculateNPscore(mol):
    """
    ---
    Ref.: Ertl, Peter, Silvio Roggo, and Ansgar Schuffenhauer.
          J Chem Inf Model, 48(1), 68-74.
    
    Based: Https://github.com/rdkit/rdkit/tree/master/Contrib/NP_Score
    
    ---
    Brief: A function to calculate the natural product-likeness score
    
    ---
    Parameters:
        >>> mol: dkit.Chem.rdchem.Mol;
    
    Return:
        float, meant NP score
    """
    filename = os.path.join(ContriDir, 'NP_Score/publicnp.model.gz')
    fscore = npscorer.pickle.load(npscorer.gzip.open(filename)) 
    return round(npscorer.scoreMol(mol,fscore=fscore),2)

    
def GetIFG(mol):
    """
    ---
    Ref.: Ertl, Peter.
          J Cheminform, 9(1), 36.
    
    Based: https://github.com/rdkit/rdkit/tree/master/Contrib/IFG
    
    ---
    Brief: A function to compute functional groups in organic molecules
    
    ---
    Parameters:
        >>> mol: dkit.Chem.rdchem.Mol;
        
    Return:
        list of namedtuple, namedtuple('IFG', ['atomIds', 'atoms', 'type'])
        e.g.:
            [IFG(atomIds=(2,), atoms='n', type='cnc'),
             IFG(atomIds=(4, 5, 6, 7), atoms='NS(=O)=O', type='cNS(c)(=O)=O'),
             IFG(atomIds=(12,), atoms='N', type='cN'),
             IFG(atomIds=(15,), atoms='n', type='cnc')]
    ---
    """
    return identify_functional_groups(mol)


def CalculateMolVolume(mol):
    """
    ---
    Equation: 
        for single atom: Vw = 4/3*pi*rw^3, the rw is the Van der Waals radius of atom
        VvdW = ∑(atom contributions)-5.92NB(Unit in Å^3), NB is the total number of bonds
        the Van der Waals radius of atom is derived from wikipedia.
        
    ---
    Parameters:
        >>> mol: dkit.Chem.rdchem.Mol;
        
    Rerurn:
        float 
    """
    from math import pi
    Radii = {'H':1.20,'C':1.70,'N':1.55,
             'O':1.52,'S':1.80,'P':1.80,
             'F':1.47,'Cl':1.75,'Br':1.85,
             'I':1.98,'Na':2.27,'Mg':1.73,
             'K':2.75,'Ca':2.31,'Ba':2.68,
             }
    mol = Chem.AddHs(mol)
    contrib = [Radii[atom.GetSymbol()] for atom in mol.GetAtoms()]
    contrib = [pi*(r**3)*4/3 for r in contrib]
    vol = sum(contrib) - 5.92*len(mol.GetBonds())
    return round(vol,2)


def CalculateMolDensity(mol):
    """
    """
    MW = CalculateMolWeight(mol)
    Vol = CalculateMolVolume(mol)
    Dens = MW/Vol
    return round(Dens,2)

def CalculateMolFCharge(mol):
    """
    """
    mol = Chem.AddHs(mol)
    FChar = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    return sum(FChar)


def _CalculateNumBond(mol,btype):
    if btype == 'SINGLE':
        return len([bond for bond in mol.GetBonds() 
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE])
    elif btype == 'DOUBLE':
        return len([bond for bond in mol.GetBonds() 
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE])
    elif btype == 'TRIPLE':
        return len([bond for bond in mol.GetBonds() 
                    if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE])


def CalculateNumSinBond(mol):
    """
    ---> nSingle
    """
    return _CalculateNumBond(mol,btype='SINGLE')


def CalculateNumDouBond(mol):
    """
    ---> nDouble
    """
    return _CalculateNumBond(mol,btype='DOUBLE')


def CalculateNumTriBond(mol,btype='TRIPLE'):
    """
    ---> nTriple
    """
    return _CalculateNumBond(mol,btype='TRIPLE')



def GetProperties(mol):
    """
    Get all properties in scopy
    """
    MW = CalculateMolWeight(mol)
    Vol = CalculateMolVolume(mol)
    Dense = CalculateMolDensity(mol)
    fChar = CalculateMolFCharge(mol)
    nBond = CalculateNumBonds(mol)
    nAtom = CalculateNumAtoms(mol)
    nCarbon = CalculateCarbonNumber(mol)
    nHet = CalculateHeteroNumber(mol)
    nRot = CalculateNumRotatableBonds(mol)
    nRig = CalculateNumRigidBonds(mol)
    nRing = CalculateNumRing(mol)
    nHev = CalculateNumHeavyAtom(mol)
    logP = CalculateLogP(mol)
    logD = CalculateLogD(mol)
    pKa = CalculatepKa(mol)
    ab = CheckAcid(mol)
    MR = CalculateMolMR(mol)
    nHD = CalculateNumHDonors(mol)
    nHA = CalculateNumHAcceptors(mol)
    nHB = CalculateNumHyBond(mol)
    AP = CalculateAromaticProportion(mol)
    logSw = CalculateLogSw(mol)
    Fsp3 = CalculateFsp3(mol)
    tPSA = CalculateTPSA(mol)
    MaxRing = CalculateMaxSizeSystemRing(mol)
    nStero = CalculateNumberStereocenters(mol)
    HetRation = CalculateHetCarbonRatio(mol)
    QED = CalculateQED(mol)
    SAscore = CalculateSAscore(mol)
    NPscore = CalculateNPscore(mol)
    nSingle = CalculateNumSinBond(mol)
    nDouble = CalculateNumDouBond(mol)
    nTriple = CalculateNumTriBond(mol)
    nC = CalculateCarbonNumber(mol)
    nB = CalculateBoronNumber(mol)
    nF = CalculateFluorinNumber(mol)
    nCl = CalculateChlorinNumber(mol)
    nBr = CalculateBromineNumber(mol)
    nI = CalculateIodineNumber(mol)
    nP = CalculatePhosphorNumber(mol)
    nS = CalculateSulfurNumber(mol)
    nO = CalculateOxygenNumber(mol)
    nN = CalculateNitrogenNumber(mol)
    
    res = namedtuple('Properties',['MW','Vol','Dense','fChar','nBond','nAtom','nCarbon','nHD','nHA','nHB',
                                   'nHet','nStero','nHev','nRot','nRig','nRing',
                                   'logP','logD','pKa','logSw','ab','MR','tPSA','AP','HetRatio',
                                   'Fsp3','MaxRing','QED','SAscore','NPscore',
                                   'nSingle','nDouble','nTriple','nC','nB','nF','nCl','nBr','nI',
                                   'nP','nS','nO','nN'])
    checkres = res(MW,Vol,Dense,fChar,nBond,nAtom,nCarbon,nHD,nHA,nHB,
                   nHet,nStero,nHev,nRot,nRig,nRing,
                   logP,logD,pKa,logSw,ab,MR,tPSA,AP,HetRation,
                   Fsp3,MaxRing,QED,SAscore,NPscore,
                   nSingle,nDouble,nTriple,nC,nB,nF,nCl,nBr,nI,
                   nP,nS,nO,nN)  
    return checkres



if __name__ =='__main__':
    
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]','CC(=O)OC1=CC=CC=C1C(=O)O']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
    smiring = ['C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2','C1CCC2C3CC(C4CCCC5CCCCC45)CCC3CCC2C1','C1CCC2(CCC3(CCCCC3)CC2)CC1']
    for index, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        print('Index:{}'.format(index))
        res = GetProperties(mol)
        print(res)
    
#    smi = 'CC(=O)OCOC(=O)C'
#    mol = Chem.MolFromSmiles(smi)
#    res = CalculateNumTriBond(mol)
#    print(res)
        
        
