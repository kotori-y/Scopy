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
        Molcular Volume >>> MolVol
        Molcular Density >>> Dense
        Formal Charge >>> fChar
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
        QED with average descriptor weights >>> QEDmean
        QED with maximal descriptor weights >>> QEDmax
        QED with using unit weights >>> QEDnone
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
        Number of double bobds >>> nDouble
        Number of triple bonds >>> nTriple
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
        Number of Charged Groups
    """
from rdkit import RDConfig
import sys,csv
sys.path.append(RDConfig.RDContribDir)
from SA_Score import sascorer
from NP_Score import npscorer
from IFG.ifg import identify_functional_groups

try:
    from ..fingerprint.fingerprints import CalculateGhoseCrippen
except:
    sys.path.append('..')
    from fingerprint.fingerprints import CalculateGhoseCrippen

try:
    from .. import ScoConfig
except:
    sys.path.append('..')
    import ScoConfig

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors, Lipinski, QED
from rdkit.Chem.Scaffolds import MurckoScaffold
from itertools import combinations
from collections import namedtuple



def CalculateMolWeight(mol):    
    """
    Calculation of molecular weight(contain hydrogen atoms)   
    ---->MW  
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    MW: float
    """
    mol = Chem.AddHs(mol)
    MW = round(sum([atom.GetMass() for atom in mol.GetAtoms()]),2)
    return MW


def CalculateNumBonds(mol):
    """
    Calculation the number of bonds where between heavy atoms       
    --->nBond
        
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nBond: int
    """
    nBond = mol.GetNumBonds()    
    return nBond


def CalculateNumAtoms(mol):
    """
    Calculation of the number of atoms in molecular 
    ---->nAtom
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nBond: int
    """  
    mol = Chem.AddHs(mol)
    return mol.GetNumAtoms()
   

def CalculateNumHetero(mol):
    """
    Calculation of Hetero counts in a molecule  
    ---->nHet
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nHet: int
    """
    i = len(
            [atom for atom in mol.GetAtoms()\
             if atom.GetAtomicNum() in [1,6]]
            )
    nHet = mol.GetNumAtoms()-i
    return nHet


def CalculateNumRotatableBonds(mol):
    """
    Calculation of the number of rotatableBonds
    ---->nRot
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nRot: int
    
    Note:
    -----------
    In some situaion Amide C-N bonds are not considered 
    because of their high rotational energy barrier
    """
    patt = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    nRot = len(mol.GetSubstructMatches(patt))    
    return nRot


def CalculateNumRigidBonds(mol):
    """
    Number of non-flexible bonds, in opposite to rotatable bonds    
    ---->nRig
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nRot: int
    """
    nBOND = CalculateNumBonds(mol)    
    flex = 0
    for bond in mol.GetBonds():
        bondtype = bond.GetBondType()
        if bondtype == Chem.rdchem.BondType.SINGLE and not bond.IsInRing():
            flex+=1 
    nRig = nBOND-flex                
    return nRig


def CalculateNumRing(mol):
    """
    Calculation of the number of ring   
    ---->nRing
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nRing: int
    """
    nRing = Chem.GetSSSR(mol)    
    return nRing


def CalculateNumHeavyAtom(mol):
    """
    Calculation of Heavy atom counts in a molecule   
    ---->nHev
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nHev: int
    """
    nHev = mol.GetNumHeavyAtoms()
    return nHev


def CalculateLogD(mol):
    """
    Calculation of molecular logD under pH=7.4
    ---->LogD
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    logD: float
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
    Calculation of molecular LogP
    ---->logP
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    logP: float
    """  
    logP = round(Descriptors.MolLogP(mol),2)  
    return logP


def CheckAcid(mol):
    """
    Judge a molecular whether is acid via SMARTS
    These SMARTS retrived from https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
    ---->ab
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    ab: str  
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
    Calculating pKa based on the ralation between logD and logP in specific pH   
    Eq.:
        |pH-pKa| = log10(10^(logP-logD)-1)
        pKa = pH - log10(10^(logP-logD)-1) for acid
        pKa = log10(10^(logP-logD)-1) - pH for base
        
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    pKa: float
    
    Note:
    -----------
    This function should be revised
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
    Cacluation of molecular refraction value based on Crippen method 
    ---->MR
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    MR: float
    """
    MR = round(Descriptors.MolMR(mol),2) 
    return MR


def CalculateNumHDonors(mol):    
    """
    Caculation of the number of Hydrogen Bond Donors
    ---->nHD
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nHD: int
    """
    nHD = Lipinski.NumHDonors(mol)    
    return nHD


def CalculateNumHAcceptors(mol):    
    """
    Caculation of the number of Hydrogen Bond Acceptors  
    ---->nHA
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nHA: int 
    """
    nHA = Lipinski.NumHAcceptors(mol)    
    return nHA


def CalculateNumHyBond(mol):
    """
    Sum of Hydrogen Bond Donnors and Acceptors   
    ---->nHB
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nHA: int 
    """
    nHD = CalculateNumHDonors(mol)
    nHA = CalculateNumHAcceptors(mol)
    nHB = nHD+nHA
    return nHB


def CalculateAromaticProportion(mol):
    """
    The proportion of heavy atoms in the molecule that are in an aromatic ring  
    ---->AP
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    AP: float 
    """
    aroma = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a]')))
    total = CalculateNumHeavyAtom(mol)
    AP = round(aroma/total,2) 
    return AP  
    

def CalculateLogSw(mol):
    """
    The logSw represents the logarithm of compounds water solubility computed by the ESOL method
    ---->logSw
    Equation: Log(Sw) = 0.16-0.638*clogP-0.0062*MWT+0.066*RB-0.74*AP
              >MWT: Molecular Weight
              >RB: Rotatable bonds
              >AP: Aromatic proportion
    
    Ref.:
    -----------
    Delaney, John S. 
    Journal of chemical information and computer sciences 44.3 (2004): 1000-1005.
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    logSw: float          
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
    Fsp3 (carbon bond saturation) is defined as the number of sp3 hybridized carbons / total carbon count.   
    ---->FSP3
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    Fsp3: float
    """
    return round(Lipinski.FractionCSP3(mol),2)
    

def CalculateTPSA(mol):
    """
    Calculation of TPSA   
    ---->TPSA
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    TPSA: float
    """
    TPSA = round(Descriptors.TPSA(mol),2)
    return TPSA
    

            
def CalculateQEDmean(mol):
    """
    Calculation QED descriptor under different weights   
    A descriptor a measure of drug-likeness based on the concept of desirability
    Here, calculating the QED descriptor using average descriptor weights.
    ---->QEDmean
    
    Ref.:
    -----------
    Bickerton, G. Richard, et al.
    Nat Chem, 4.2 (2012): 90.
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    qed: float
    """    
    QEDmean = QED.weights_mean(mol)        
    return round(QEDmean,2) 


def CalculateQEDmax(mol):
    """
    Calculation QED descriptor under different weights   
    A descriptor a measure of drug-likeness based on the concept of desirability
    Here, calculating the QED descriptor using maximal descriptor weights.
    ---->QEDmax
    
    Ref.:
    -----------
    Bickerton, G. Richard, et al.
    Nat Chem, 4.2 (2012): 90.
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    qed: float
    """    
    QEDmax = QED.weights_max(mol)        
    return round(QEDmax,2)     


def CalculateQEDnone(mol):
    """
    Calculation QED descriptor under different weights   
    A descriptor a measure of drug-likeness based on the concept of desirability
    Here, calculating the QED descriptor using unit weights.
    ---->QEDnone
    
    Ref.:
    -----------
    Bickerton, G. Richard, et al.
    Nat Chem, 4.2 (2012): 90.
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    qed: float
    """    
    QEDnone = QED.weights_none(mol)        
    return round(QEDnone,2)


def CalculateMaxSizeSystemRing(mol):
    """
    Number of atoms involved in the biggest system ring  
    ----> maxring
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    MaxRing: int
    """
    #0.Get the scaffold
    core = MurckoScaffold.GetScaffoldForMol(mol)
    fw = MurckoScaffold.MakeScaffoldGeneric(core)
    #1.Obtaining which atoms consist of rings
    MaxRing = 0
    ri = fw.GetRingInfo()
    atoms = list(ri.AtomRings())    
    length = len(atoms)    
    if length == 0:
        pass
    else:
        rw = Chem.RWMol(fw)        
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
            if len(item) > MaxRing:
                MaxRing = len(item)   
    return MaxRing
    

def CalculateNumberStereocenters(mol):
    """
    the number of stereo centers
    ---->nStero
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nStero: int
    """
    return Chem.CalcNumAtomStereoCenters(mol)    


def _CalculateElementNumber(mol,AtomicNumber=6):
    """
    **Internal used only**
    Calculation of specific type of atom number in a molecule
    
    Calculation of element counts with atomic number equal to n in a molecule
    """
    return len(
            [atom for atom in mol.GetAtoms()\
                if atom.GetAtomicNum() == AtomicNumber]
            )


def CalculateCarbonNumber(mol):
    """
    Calculation of Carbon number in a molecule    
    ---->nC
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nC: int
    """
    return _CalculateElementNumber(mol,AtomicNumber=6)


def CalculateBoronNumber(mol):
    """
    Calculation of Boron counts in a molecule  
    ---->nB
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nB: int
    """       
    return _CalculateElementNumber(mol,AtomicNumber=5)


def CalculateFluorinNumber(mol):
    """
    Calculation of Fluorin counts in a molecule  
    ---->nF
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nF: int
    """         
    return _CalculateElementNumber(mol,AtomicNumber=10)


def CalculateChlorinNumber(mol):
    """
    Calculation of Chlorin counts in a molecule
    ---->nCl
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nCl: int
    """
    return _CalculateElementNumber(mol,AtomicNumber=17)


def CalculateBromineNumber(mol):
    """
    Calculation of Bromine counts in a molecule  
    ---->nBr
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nBr: int
    """
    return _CalculateElementNumber(mol,AtomicNumber=35)


def CalculateIodineNumber(mol):
    """
    Calculation of Iodine counts in a molecule 
    ---->nI
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nI: int
    """
    return _CalculateElementNumber(mol,AtomicNumber=53)


def CalculatePhosphorNumber(mol):
    """
    Calcualtion of Phosphor number in a molecule
    ---->nP
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nP: int
    """
    return _CalculateElementNumber(mol,AtomicNumber=15)


def CalculateSulfurNumber(mol):
    """
    Calculation of Sulfur counts in a molecule  
    ---->nS
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nCa: int
    """
    return _CalculateElementNumber(mol,AtomicNumber=16)


def CalculateOxygenNumber(mol):
    """
    Calculation of Oxygen counts in a molecule    
    ---->nO
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nO: int
    """
    return _CalculateElementNumber(mol,AtomicNumber=8)
        

def CalculateNitrogenNumber(mol):
    """
    Calculation of Nitrogen counts in a molecule
    ---->nnitro
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nN: int
    """
    return _CalculateElementNumber(mol,AtomicNumber=7)


def CalculateNumberChargedGroups(mol):
    """
    Number of Charged Groups 
    ---->nChar
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nChar: int
    """
    pass


def CalculateHetCarbonRatio(mol):
    """
    The ratio between the number of non carbon atoms and the number of carbon atoms.
    ---->HetRatio
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    HetRatio: float
    """
    Total = CalculateNumHeavyAtom(mol)
    nCarb = CalculateCarbonNumber(mol)
    het = Total-nCarb
    try:
        HetRatio = round(het/nCarb,2)
    except:
        HetRatio = 'Inf'
    return HetRatio   
    

def CalculateSAscore(mol):
    """
    A function to estimate ease of synthesis (synthetic accessibility) of drug-like molecules
    ---->SAscore
    
    Ref.:
    -----------
    Ertl, Peter, and Ansgar Schuffenhauer.
    J Cheminform, 1(1), 8.      
    Based: https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score   
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    SAscore: float
    """
    return round(sascorer.calculateScore(mol),2)


def CalculateNPscore(mol):
    import os
    """
    A function to calculate the natural product-likeness score
    ---->NPscore
    
    Ref.:
    -----------
    Ertl, Peter, Silvio Roggo, and Ansgar Schuffenhauer.
    J Chem Inf Model, 48(1), 68-74.
    Based: Https://github.com/rdkit/rdkit/tree/master/Contrib/NP_Score
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    NPscore: float
    """
    ContriDir = RDConfig.RDContribDir
    filename = os.path.join(ContriDir, 'NP_Score/publicnp.model.gz')
    fscore = npscorer.pickle.load(npscorer.gzip.open(filename)) 
    return round(npscorer.scoreMol(mol,fscore=fscore),2)

    
def GetIFG(mol):
    """
    A function to compute functional groups in organic molecules
    
    Ref.: 
    -----------
    Ertl, Peter.
    J Cheminform, 9(1), 36.
    Based: https://github.com/rdkit/rdkit/tree/master/Contrib/IFG
    
    Parameters:
    -----------
    mol: dkit.Chem.rdchem.Mol;
        
    Return:
    -----------
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
    Calculation of Van der Waals Volume of molecule
    ---->MolVol
    
    Equation: 
    -----------
    for single atom: Vw = 4/3*pi*rw^3, the rw is the Van der Waals radius of atom
    VvdW = ∑(atom contributions)-5.92NB(Unit in Å^3), NB is the total number of bonds
    the Van der Waals radius of atom is derived from wikipedia.
        
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    NPscore: float
    """
    from math import pi
    Radii = {'H':1.20,'C':1.70,'N':1.55,
             'O':1.52,'S':1.80,'P':1.80,
             'F':1.47,'Cl':1.75,'Br':1.85,
             'I':1.98,'Na':2.27,'Mg':1.73,
             'K':2.75,'Ca':2.31,'Ba':2.68,
             'He':140,'Li':182,'Be':153,
             'B':192,'Ne':154,'Al':184,
             'Si':210,'Ar':188,'Ni':163,
             'Cu':140,'Zn':139,'Ga':187,
             'Ge':211,'As':185,'Se':190,
             'Kr':202,'Rb':303,'Sr':249,
             'Pd':163,'Ag':172,'Cd':158,
             'In':193,'Sn':217,'Sb':206,
             'Te':206,'Xe':216,'Cs':343,
             'Pt':175,'Au':166,'U':186,
             'Hg':155,'Tl':196,'Pb':202,
             'Bi':207,'Po':197,'At':202,
             'Rn':220,'Fr':348,'Ra':283}
    mol = Chem.AddHs(mol)
    contrib = []
    for atom in mol.GetAtoms():
        try:
            contrib.append(Radii[atom.GetSymbol()])
        except:
            pass
    # contrib = [Radii[atom.GetSymbol()] for atom in mol.GetAtoms()]
    contrib = [pi*(r**3)*4/3 for r in contrib]
    vol = sum(contrib) - 5.92*len(mol.GetBonds())
    return round(vol,2)


def CalculateMolDensity(mol):
    """
    Calculation of Density of molecule
    ---->Dense
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    Dense: float
    """
    MW = CalculateMolWeight(mol)
    Vol = CalculateMolVolume(mol)
    Dens = MW/Vol
    return round(Dens,2)

def CalculateMolFCharge(mol):
    """
    Calculation of formal charge of molecule
    ---->fChar
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    fChar: float
    """
    mol = Chem.AddHs(mol)
    FChar = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    return sum(FChar)


def _CalculateNumBond(mol,btype):
    """
    **Internal used only**
    Calculation of specific type of bond number in a molecule
    """
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
    Calculation of single bond number of molecule
    ---> nSingle
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nSingle: int
    """
    return _CalculateNumBond(mol,btype='SINGLE')


def CalculateNumDouBond(mol):
    """
    Calculation of double bond number of molecule
    --->nDouble
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nDouble: int
    """
    return _CalculateNumBond(mol,btype='DOUBLE')


def CalculateNumTriBond(mol,btype='TRIPLE'):
    """
    Calculation of triple bond number of molecule
    ---> nTriple
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    
    Return:
    -----------
    nTriple: int
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
    nHet = CalculateNumHetero(mol)
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
    HetRatio = CalculateHetCarbonRatio(mol)
    QEDmean = CalculateQEDmean(mol)
    QEDmax = CalculateQEDmax(mol)
    QEDnone = CalculateQEDnone(mol)
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
                                   'Fsp3','MaxRing','QEDmean','QEDmax','QEDnone','SAscore','NPscore',
                                   'nSingle','nDouble','nTriple','nC','nB','nF','nCl','nBr','nI',
                                   'nP','nS','nO','nN'])
    checkres = res(MW,Vol,Dense,fChar,nBond,nAtom,nCarbon,nHD,nHA,nHB,
                   nHet,nStero,nHev,nRot,nRig,nRing,
                   logP,logD,pKa,logSw,ab,MR,tPSA,AP,HetRatio,
                   Fsp3,MaxRing,QEDmean,QEDmax,QEDnone,SAscore,NPscore,
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
        res = GetIFG(mol)
        print(res)
    
#    smi = 'CC(=O)OCOC(=O)C'
#    mol = Chem.MolFromSmiles(smi)
#    res = CalculateNumTriBond(mol)
#    print(res)
        
        
