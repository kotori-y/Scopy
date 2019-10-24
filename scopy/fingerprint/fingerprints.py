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
try:
    from ghosecrippen import GCfp
    from estate import EStateFP
    from morgan import Morgan
    from daylight import Daylight
    from efg import EFG
    from pubchem import calcPubChemFingerAll
except:
    from .ghosecrippen import GCfp
    from .estate import EStateFP
    from .morgan import Morgan
    from .daylight import Daylight
    from .efg import EFG
    from .pubchem import calcPubChemFingerAll

from rdkit.Chem.rdMolDescriptors import GetMACCSKeysFingerprint
import numpy as np



def CalculateEFG(mols,counter=True):
    """
    Brief:
    -----------
    classification system termed “extended functional groups” (EFG),
    which are an extension of a set previously used by the CheckMol software, 
    that covers in addition heterocyclic compound classes and periodic table groups. 
    583 bits
    
    Ref.:
    -----------
    Salmina, Elena, Norbert Haider, and Igor Tetko.
    Molecules, 21(1), 1.
    
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    counter: bool(optional, default:True)
        If set to True, the fingerprint will presented in the format of counter(not only 1 and 0)
        else, would be binary
        
    Return:
    -----------
    fps: numpy.ndarray
    """
    fps = EFG(counter=counter)
    fps = list(map(lambda mol: fps.CalculateEFG(mol), mols))
    fps = np.array(fps)
    return fps


def CalculateGhoseCrippen(mols,counter=True):
    """
    Brief:
    -----------
    Atom-based calculation of LogP and MR using Crippen’s approach, it's 110 bits
    
    Ref.:
    -----------
    Wildman, Scott A., and Gordon M. Crippen.
    JCICS, 39(5), 868-873.
    
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    counter: bool(optional, default:True)
        If set to True, the fingerprint will presented in the format of counter(not only 1 and 0)
        else, would be binary
        
    Return:
    -----------
    fps: numpy.ndarray
    """
    fps = GCfp(counter)
    fps = list(map(lambda mol: fps.CalculateGCfp(mol), mols))
    fps = np.array(fps)
    return fps
    
            
def CalculateEState(mols,val=True):
    """
    Brief:
    -----------
    Atom-based calculation of LogP and MR using Crippen’s approach, 79 bits
    
    Ref.:
    -----------
    Kier, Lemont Burwell, and Lowell H. Hall.
    Academic, 1999
    
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    val: bool(optional, default:True)
        If set to True, the fingerprint will presented in the format of vaule of estate
        else, would be binary
        
    Return:
    -----------
    fps: numpy.ndarray
    """
    fps = EStateFP(val)
    fps = list(map(lambda mol: list(fps.CalculateEState(mol)), mols))
    fps = np.array(fps)
    return fps


def CalculateMACCS(mols):
    """
    Brief:
    -----------
    There is a SMARTS-based implementation of the 166 public MACCS keys. 166 bits
    
    Ref.:
    -----------
    Using the 166 public keys implemented as SMARTS
    
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
        
    Return:
    -----------
    fps: numpy.ndarray
    """
    fps = list(map(lambda mol: np.array(GetMACCSKeysFingerprint(mol)), mols))
    fps = np.array(fps)
    return fps


def CalculateECFP(mols, radius=2, nBits=1024):
    """
    Brief:
    -----------
    This family of fingerprints, better known as circular fingerprints, 
    is built by applying the Morgan algorithm to a set of user-supplied atom invariants.
    2^n bits
    
    Ref.:
    -----------
    Rogers, David, and Mathew Hahn.
    J Chem Inf Model, 50(5), 742-754.
    
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    radius: int(optional, defailt: 2)
        the radius of circle
    nBits: int(optional, default:1024)
        number of bits in the fingerprint
        
    Return:
    -----------
    fps: numpy.ndarray
    """
    fps = Morgan(radius=radius, nBits=nBits)
    fps = list(map(lambda mol: fps.CalculateECFP(mol), mols))
    fps = np.array(fps)
    return fps


def CalculateFCFP(mols, radius=2, nBits=1024):
    """
    Brief:
    -----------
    This family of fingerprints, better known as circular fingerprints, 
    is built by applying the Morgan algorithm to a set of user-supplied atom invariants.
    2^n bits
    
    Ref.:
    -----------
    Rogers, David, and Mathew Hahn.
    J Chem Inf Model, 50(5), 742-754.
    
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    radius: int(optional, defailt: 2)
        the radius of circle
    nBits: int(optional, default:1024)
        number of bits in the fingerprint
        
    Return:
    -----------
    fps: numpy.ndarray
    """
    fps = Morgan(radius=radius, nBits=nBits)
    fps = list(map(lambda mol: fps.CalculateECFP(mol), mols))
    fps = np.array(fps)
    return fps


def CalculateDaylight(mols, minPath=1, maxPath=7, nBits=2048):
    """
    Brief:
    -----------
    a Daylight-like fingerprint based on hashing molecular subgraphs
    
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    minPath: int(optional, default:1)
        minimum number of bonds to include in the subgraphs
    maxPath: int(optional, default:7)
        maximum number of bonds to include in the subgraphs
    nBits: int(optional, default:2048)
        number of bits in the fingerprint
        
    Return:
    -----------
    fps: numpy.ndarray
    """
    fps = Daylight(minPath=minPath, maxPath=maxPath, nBits=nBits)
    fps = list(map(lambda mol: fps.CalculateDaylight(mol), mols))
    fps = np.array(fps)
    return fps

                                           
def CalculatePubChemFingerprint(mols):
    """
    Parameters:
    -----------
    mols: Iterable object, each element is a rdkit.Chem.rdchem.Mol
        The molecule(s) to be scanned
    
    Return:
    -----------
    fps: numpy.ndarray   
    """
    fps = list(map(lambda mol: calcPubChemFingerAll(mol), mols))
    fps = np.array(fps)
    return fps
 

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
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    fps = CalculateGhoseCrippen(mols)
    print(fps)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

