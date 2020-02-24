# -*- coding: utf-8 -*-

#Created on Tue Jun 25 21:59:42 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me



from itertools import combinations
import sys,csv,os

from rdkit import RDConfig
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors, Lipinski, QED
from rdkit.Chem.Scaffolds import MurckoScaffold
sys.path.append(RDConfig.RDContribDir)

from SA_Score import sascorer
from NP_Score import npscorer
from IFG.ifg import identify_functional_groups

try:
    from ..fingerprint import fingerprints
except:
    sys.path.append(__file__+'//..')
    from fingerprint import fingerprints

try:
    from .. import ScoConfig
except:
    sys.path.append('..')
    import ScoConfig



ContriDir = RDConfig.RDContribDir
filename = os.path.join(ContriDir, 'NP_Score/publicnp.model.gz')
fscore = npscorer.pickle.load(npscorer.gzip.open(filename)) 


def CalculateMolWeight(mol):    
    """
    Calculation of molecular weight(contain hydrogen atoms)   
    --->MW  

    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the weight of molecule(contain hydrogen atoms)
    :rtype: float
    
    """
    MW = Descriptors.ExactMolWt(mol)
    return round(MW,2)


def CalculateNumBonds(mol):
    """
    Calculation the number of bonds where between heavy atoms       
    --->nBond
        
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of bonds where between heavy atoms
    :rtype: int
    
    """
    nBond = mol.GetNumBonds()    
    return nBond


def CalculateNumAtoms(mol):
    """
    Calculation of the number of atoms in molecular(contain hydrogen atoms) 
    --->nAtom
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of atoms in molecular(contain hydrogen atoms)
    :rtype: int
    
    """  
    mol = Chem.AddHs(mol)
    return mol.GetNumAtoms()
   

def CalculateNumHetero(mol):
    """
    Calculation of the number of heteroatom in a molecule  
    --->nHet
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of heteroatom in a molecule  
    :rtype: int
    
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
    --->nRot
    
    Note:
        In some situaion Amide C-N bonds are not considered 
        because of their high rotational energy barrier
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of rotatableBond  
    :rtype: int
    
    
    """
    patt = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    nRot = len(mol.GetSubstructMatches(patt))    
    return nRot


def CalculateNumRigidBonds(mol):
    """
    Number of non-flexible bonds, in opposite to rotatable bonds    
    --->nRig
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of non-flexible bonds 
    :rtype: int
    
    """
    nBOND = CalculateNumBonds(mol)    
    flex = 0
    for bond in mol.GetBonds():
        bondtype = bond.GetBondType()
        if bondtype == Chem.rdchem.BondType.SINGLE and not bond.IsInRing():
            flex+=1 
    nRig = nBOND-flex                
    return nRig


def CalculateFlexibility(mol):
    """
    The flexibility (ration between rotatable and rigid bonds)
    --->Flex
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of ring   
    :rtype: float
    
    """
    nRot = CalculateNumRotatableBonds(mol)
    nRig = CalculateNumRigidBonds(mol)
    return round(nRot/nRig,2) if nRig else 'Inf'
   
    
def CalculateNumRing(mol):
    """
    Calculation of the number of ring   
    --->nRing
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of ring   
    :rtype: int
    
    """
    nRing = Chem.GetSSSR(mol)    
    return nRing


def CalculateNumHeavyAtom(mol):
    """
    Calculation of Heavy atom counts in a molecule   
    --->nHev
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of heavy atom counts in a molecule  
    :rtype: int
    
    """
    nHev = mol.GetNumHeavyAtoms()
    return nHev


def CalculateLogD(mol):
    """
    Calculation of molecular logD under pH=7.4
    --->LogD
    
    Note:
        We have built a liner model with DNN to predict logD7.4.
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: molecular logD under pH=7.4 
    :rtype: float
    
    """
    intercept = 0.5748907159915493
    
    fps = fingerprints.CalculateGhoseCrippen([mol]).flatten()
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
    --->logP
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: molecular logP 
    :rtype: float
    
    """  
    logP = round(Descriptors.MolLogP(mol),2)  
    return logP


def CheckAcid(mol):
    """
    Judge a molecular whether is acid via SMARTS.
    These SMARTS retrived from https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
    --->ab
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: classification to acid or base
    :rtype: str
    
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
    Calculating pKa based on the ralation between logD and logP in specific pH.
    --->pKa
    
    Eq.:
        abs(pH-pKa) = log10(10^(logP-logD)-1)
        pKa = pH - log10(10^(logP-logD)-1) for acid
        pKa = log10(10^(logP-logD)-1) - pH for base
        
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: molecular pKa
    :rtype: float
    
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
    --->MR
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: molecular refraction value based on Crippen method
    :rtype: float
    
    """
    MR = round(Descriptors.MolMR(mol),2) 
    return MR


def CalculateNumHDonors(mol):    
    """
    Caculation of the number of Hydrogen Bond Donors
    --->nHD
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of Hydrogen Bond Donors
    :rtype: int
    
    """
    nHD = Lipinski.NumHDonors(mol)    
    return nHD


def CalculateNumHAcceptors(mol):    
    """
    Caculation of the number of Hydrogen Bond Acceptors  
    --->nHA
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of Hydrogen Bond Acceptors
    :rtype: int
    
    """
    nHA = Lipinski.NumHAcceptors(mol)    
    return nHA


def CalculateNumHyBond(mol):
    """
    Sum of Hydrogen Bond Donnors and Acceptors   
    --->nHB
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: sum of Hydrogen Bond Donnors and Acceptors
    :rtype: int
    
    """
    nHD = CalculateNumHDonors(mol)
    nHA = CalculateNumHAcceptors(mol)
    nHB = nHD+nHA
    return nHB


def CalculateNumAromaAtom(mol):
    """
    Calculation of aromatic atom counts in a molecule
    --->nAAtom
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the aromatic atom counts in a molecule
    :rtype: int
    
    """
    aroma = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a]')))
    return aroma


def CalculateNumAromaRing(mol):
    n = 0
    aatom = mol.GetSubstructMatches(Chem.MolFromSmarts('[a]'))
    aatom = sum(aatom,())
    if aatom:
        ringinfo = mol.GetRingInfo()
        for info in ringinfo.AtomRings():
            n += 1 if set(info) == set(info)&set(aatom) else 0
    else:
        pass
    return n


def CalculateAromaticProportion(mol):
    """
    The proportion of heavy atoms in the molecule that are in an aromatic ring  
    --->AP
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the proportion of heavy atoms in the molecule that are in an aromatic ring  
    :rtype: float
    
    """
    aroma = CalculateNumAromaAtom(mol)
    total = CalculateNumHeavyAtom(mol)
    AP = round(aroma/total,2) 
    return AP  
    

def CalculateLogSw(mol):
    """
    The logSw represents the logarithm of compounds water solubility computed by the ESOL method
    --->logSw
    
    Equation: 
        Log(Sw) = 0.16-0.638*clogP-0.0062*MWT+0.066*RB-0.74*AP
        where, MWT: Molecular Weight; RB: Rotatable bonds; AP: Aromatic proportion
    
    Reference:
        (1) `Delaney, John S (2004)`_. 
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the molecular logSw
    :rtype: float
    
    .. _Delaney, John S (2004):
        https://pubs.acs.org/doi/abs/10.1021/ci034243x
        
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
    --->FSP3
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the carbon bond saturation
    :rtype: float
    
    """
    return round(Lipinski.FractionCSP3(mol),2)
    

def CalculateTPSA(mol):
    """
    Calculation of TPSA   
    --->TPSA
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: TPSA
    :rtype: float
    
    """
    TPSA = round(Descriptors.TPSA(mol),2)
    return TPSA
    
           
def CalculateQEDmean(mol):
    """
    Calculation QED descriptor under different weights 
    A descriptor a measure of drug-likeness based on the concept of desirability
    Here, calculating the QED descriptor using average descriptor weights.
    --->QEDmean
    
    Reference:
        (1) `Bickerton, G. Richard (2012)`_.
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: QED descriptor using average descriptor weights
    :rtype: float
    
    .. _Bickerton, G. Richard (2012):
        https://www.nature.com/nchem/journal/v4/n2/abs/nchem.1243.html
        
    """    
    QEDmean = QED.weights_mean(mol)        
    return round(QEDmean,2) 


def CalculateQEDmax(mol):
    """
    Calculation QED descriptor under different weights   
    A descriptor a measure of drug-likeness based on the concept of desirability
    Here, calculating the QED descriptor using maximal descriptor weights.
    --->QEDmax
    
    Reference:
        (1) `Bickerton, G. Richard (2012)`_.
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: QED descriptor using maximal descriptor weights
    :rtype: float
    
    .. _Bickerton, G. Richard (2012):
        https://www.nature.com/nchem/journal/v4/n2/abs/nchem.1243.html
        
    """    
    QEDmax = QED.weights_max(mol)        
    return round(QEDmax,2)     


def CalculateQEDnone(mol):
    """
    Calculation QED descriptor under different weights   
    A descriptor a measure of drug-likeness based on the concept of desirability
    Here, calculating the QED descriptor using unit weights.
    --->QEDnone
    
    Reference:
        (1) `Bickerton, G. Richard (2012)`_.
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: QED descriptor using unit weights
    :rtype: float
    
    .. _Bickerton, G. Richard (2012):
        https://www.nature.com/nchem/journal/v4/n2/abs/nchem.1243.html
        
    """    
    QEDnone = QED.weights_none(mol)        
    return round(QEDnone,2)


def CalculateMaxSizeSystemRing(mol):
    """
    Number of atoms involved in the biggest system ring  
    ---> maxring
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: number of atoms involved in the biggest system ring
    :rtype: int
    
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
    

def CalculateNumStereocenters(mol):
    """
    the number of stereo centers
    --->nStereo
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of stereo centers
    :rtype: int
    
    """
    return Chem.CalcNumAtomStereoCenters(mol)    


def _CalculateNumElement(mol,AtomicNumber=6):
    """
    **Internal used only**
    Calculation of specific type of atom number in a molecule
    
    Calculation of element counts with atomic number equal to n in a molecule
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :param AtomicNumber: the AtomicNumber of atom to be counted, defaults to 6
    :type AtomicNumber: int, optional
    :return: the number of stereo centers
    :rtype: int
    
    """
    return len(
            [atom for atom in mol.GetAtoms()\
                if atom.GetAtomicNum() == AtomicNumber]
            )


def CalculateNumCarbon(mol):
    """
    Calculation of Carbon number in a molecule    
    --->nC
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of carbon atoms
    :rtype: int
    
    """
    return _CalculateNumElement(mol,AtomicNumber=6)


def CalculateNumBoron(mol):
    """
    Calculation of Boron counts in a molecule  
    --->nB
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of boron atoms
    :rtype: int
    
    """       
    return _CalculateNumElement(mol,AtomicNumber=5)


def CalculateNumFluorin(mol):
    """
    Calculation of Fluorin counts in a molecule  
    --->nF
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of fluori atoms
    :rtype: int
    
    """         
    return _CalculateNumElement(mol,AtomicNumber=10)


def CalculateNumChlorin(mol):
    """
    Calculation of Chlorin counts in a molecule
    --->nCl
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of chlorin atoms
    :rtype: int
    
    """
    return _CalculateNumElement(mol,AtomicNumber=17)


def CalculateNumBromine(mol):
    """
    Calculation of Bromine counts in a molecule  
    --->nBr
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of bromine atoms
    :rtype: int
    
    """
    return _CalculateNumElement(mol,AtomicNumber=35)


def CalculateNumIodine(mol):
    """
    Calculation of Iodine counts in a molecule 
    --->nI
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of bromine atoms
    :rtype: int
    
    """
    return _CalculateNumElement(mol,AtomicNumber=53)


def CalculateNumPhosphor(mol):
    """
    Calcualtion of Phosphor number in a molecule
    --->nP
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of phosphor atoms
    :rtype: int
    
    """
    return _CalculateNumElement(mol,AtomicNumber=15)


def CalculateNumSulfur(mol):
    """
    Calculation of Sulfur counts in a molecule  
    --->nS
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of sulfur atoms
    :rtype: int
    
    """
    return _CalculateNumElement(mol,AtomicNumber=16)


def CalculateNumOxygen(mol):
    """
    Calculation of Oxygen counts in a molecule    
    --->nO
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of oxygen atoms
    :rtype: int
    
    """
    return _CalculateNumElement(mol,AtomicNumber=8)
        

def CalculateNumNitrogen(mol):
    """
    Calculation of Nitrogen counts in a molecule
    --->nN
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of nitrogen atoms
    :rtype: int
    
    """
    return _CalculateNumElement(mol,AtomicNumber=7)


# def CalculateNumberChargedGroups(mol):
#     """
#     Number of Charged Groups 
#     --->nChar
    
#     :param mol: molecular
#     :type mol: rdkit.Chem.rdchem.Mol
#     :return: the number of charged group
#     :rtype: int
    
#     """
#     pass


def CalculateHetCarbonRatio(mol):
    """
    The ratio between the number of non carbon atoms and the number of carbon atoms.
    --->HetRatio
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the ratio between the number of non carbon atoms and the number of carbon atoms
    :rtype: float
    
    """
    nHet = CalculateNumHetero(mol)
    nCarb = CalculateNumCarbon(mol)
    return round(nHet/nCarb,2) if nCarb else 'Inf'
    

def CalculateSAscore(mol):
    """
    A function to estimate ease of synthesis (synthetic accessibility) of drug-like molecules
    --->SAscore
    
    Reference:
        (1) `Ertl Peter (2009)`_.
        
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: ease of synthesis
    :rtype: float
    
    .. _Ertl Peter (2009):
        https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8
        
    """
    return round(sascorer.calculateScore(mol),2)


def CalculateNPscore(mol):
    """
    A function to calculate the natural product-likeness score
    --->NPscore
    
    Reference:
        (1) `Ertl Peter (2008)`_.
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: product-likeness score
    :rtype: float
    
    .. _Ertl Peter (2008):
        https://pubs.acs.org/doi/abs/10.1021/ci700286x
    
    """
    return round(npscorer.scoreMol(mol,fscore=fscore),2)

    
def GetIFG(mol):
    """
    A function to compute functional groups in organic molecules
    --->IFG
    
    Reference:
        (1) `Ertl Peter (2017)`_.
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: list of namedtuple, namedtuple('IFG', ['atomIds', 'atoms', 'type'])
    :rtype: list
    
    .. _Ertl Peter (2017):
        https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0225-z
        
    """
    return [fg._asdict() for fg in identify_functional_groups(mol)]


def CalculateMolVolume(mol):
    """
    Calculation of Van der Waals Volume of molecule
    --->MV
    
    Equation: 
        for single atom: Vw = 4/3*pi*rw^3, the rw is the Van der Waals radius of atom
        VvdW = ∑(atom contributions)-5.92NB(Unit in Å^3), NB is the total number of bonds
        the Van der Waals radius of atom is derived from wikipedia.
        
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: Van der Waals Volume of molecule
    :rtype: float
    
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
    Calculation of density of molecule
    --->Dense
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: density of molecule
    :rtype: float
    
    """
    MW = CalculateMolWeight(mol)
    Vol = CalculateMolVolume(mol)
    return round(MW/Vol, 2) if Vol else 'Inf'


def CalculateMolFCharge(mol):
    """
    Calculation of formal charge of molecule
    --->fChar
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: formal charge of molecule
    :rtype: float
    
    """
    mol = Chem.AddHs(mol)
    FChar = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    return sum(FChar)


def _CalculateNumBond(mol,btype):
    """
    **Internal used only**
    Calculation of specific type of bond number in a molecule
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :param btype: the type of bond to be counted
    :type btype: str
    :return: the number of specific bond
    :rtype: int
    
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
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of single bond
    :rtype: int
    
    """
    return _CalculateNumBond(mol,btype='SINGLE')


def CalculateNumDouBond(mol):
    """
    Calculation of double bond number of molecule
    --->nDouble
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of double bond
    :rtype: int
    
    """
    return _CalculateNumBond(mol,btype='DOUBLE')


def CalculateNumTriBond(mol,btype='TRIPLE'):
    """
    Calculation of triple bond number of molecule
    ---> nTriple
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: the number of triple bond
    :rtype: int
    
    """
    return _CalculateNumBond(mol,btype='TRIPLE')


def GetProperties(mol, 
                  items=['MW','Vol','Dense','fChar','nBond','nAtom','nHD','nHA','nHB',
                         'nHet','nStero','nHev','nRot','nRig','Flex','nRing','logP',
                         'logD','pKa','logSw','ab','MR','TPSA','AP','HetRatio','Fsp3',
                         'MaxRing','QEDmean','QEDmax','QEDnone','SAscore','NPscore',
                         'nSingle','nDouble','nTriple','nC','nB','nF','nCl','nBr','nI',
                         'nP','nS','nO','nN']
        ):
    """
    Get all properties in scopy
    """
    funcl = {'MW': 'CalculateMolWeight(mol)',
    'Vol': 'CalculateMolVolume(mol)',
    'Dense': 'CalculateMolDensity(mol)',
    'fChar': 'CalculateMolFCharge(mol)',
    'nBond': 'CalculateNumBonds(mol)',
    'nAtom': 'CalculateNumAtoms(mol)',
    'nHet': 'CalculateNumHetero(mol)',
    'nRot': 'CalculateNumRotatableBonds(mol)',
    'nRig': 'CalculateNumRigidBonds(mol)',
    'Flex': 'CalculateFlexibility(mol)',
    'nRing': 'CalculateNumRing(mol)',
    'nHev': 'CalculateNumHeavyAtom(mol)',
    'logP': 'CalculateLogP(mol)',
    'logD': 'CalculateLogD(mol)',
    'pKa': 'CalculatepKa(mol)',
    'ab': 'CheckAcid(mol)',
    'MR': 'CalculateMolMR(mol)',
    'nHD': 'CalculateNumHDonors(mol)',
    'nHA': 'CalculateNumHAcceptors(mol)',
    'nHB': 'CalculateNumHyBond(mol)',
    'AP': 'CalculateAromaticProportion(mol)',
    'logSw': 'CalculateLogSw(mol)',
    'Fsp3': 'CalculateFsp3(mol)',
    'TPSA': 'CalculateTPSA(mol)',
    'MaxRing': 'CalculateMaxSizeSystemRing(mol)',
    'nStero': 'CalculateNumStereocenters(mol)',
    'HetRatio': 'CalculateHetCarbonRatio(mol)',
    'QEDmean': 'CalculateQEDmean(mol)',
    'QEDmax': 'CalculateQEDmax(mol)',
    'QEDnone': 'CalculateQEDnone(mol)',
    'SAscore': 'CalculateSAscore(mol)',
    'NPscore': 'CalculateNPscore(mol)',
    'nSingle': 'CalculateNumSinBond(mol)',
    'nDouble': 'CalculateNumDouBond(mol)',
    'nTriple': 'CalculateNumTriBond(mol)',
    'nC': 'CalculateNumCarbon(mol)',
    'nB': 'CalculateNumBoron(mol)',
    'nF': 'CalculateNumFluorin(mol)',
    'nCl': 'CalculateNumChlorin(mol)',
    'nBr': 'CalculateNumBromine(mol)',
    'nI': 'CalculateNumIodine(mol)',
    'nP': 'CalculateNumPhosphor(mol)',
    'nS': 'CalculateNumSulfur(mol)',
    'nO': 'CalculateNumOxygen(mol)',
    'nN': 'CalculateNumNitrogen(mol)'}
    
    vals = []
    for item in items:
        val = eval(funcl[item])
        vals.append(val)
        
    return dict(zip(items, vals))




if __name__ =='__main__':
    
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
    
    for index, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        print('Index:{}'.format(index))
        res = GetProperties(mol=mol)
        print(res)
    

        
        
