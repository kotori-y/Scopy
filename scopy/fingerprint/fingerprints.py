# -*- coding: utf-8 -*-

#Created on Tue Jun 25 21:59:42 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me


from rdkit.Chem import AllChem as Chem
try:
    from .ghosecrippen import GCfp
    from .estate import EStateFP
    from .morgan import Morgan
    from .daylight import Daylight
    from .efg import EFG
    from .pubchem import PubChem
    from .maccs import MACCS
    from .ifg import IFG
except:
    import sys
    sys.path.append('.')
    from ghosecrippen import GCfp
    from estate import EStateFP
    from morgan import Morgan
    from daylight import Daylight
    from efg import EFG
    from pubchem import PubChem
    from maccs import MACCS
    from ifg import IFG
import numpy as np
from multiprocessing import Pool


def CalculateEFG(mols, useCount=True, n_jobs=1):
    """classification system termed “extended functional groups” (EFG),
    which are an extension of a set previously used by the CheckMol software, 
    that covers in addition heterocyclic compound classes and periodic table groups. 
    583 bits
    
    Reference:
        (1) `Salmina, Elena, Norbert Haider and Igor Tetko (2016)`_.
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param useCount: If set to True, the fingerprint will presented in the format of counter(not only 1 and 0) else, would be binary, defaults to True
    :type useCounter: bool, optional
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
    
    .. _Salmina, Elena, Norbert Haider and Igor Tetko (2016):
        https://www.mdpi.com/1420-3049/21/1/1
        
    """
    n_jobs = n_jobs if n_jobs >= 1 else None
    fps = EFG(useCount=useCount)
    pool = Pool(n_jobs)
    fps = pool.map_async(fps.CalculateEFG, mols).get()
    pool.close()
    pool.join()
    fps = np.array(fps)
    return fps


def CalculateGhoseCrippen(mols, useCount=True, n_jobs=1):
    """Atom-based calculation of LogP and MR using Crippen’s approach. 
    110 bits
   
    Reference:
        (1) `Wildman, Scott A., and Gordon M. Crippen (1999)`_.
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param useCount: If set to True, the fingerprint will presented in the format of counter(not only 1 and 0) else, would be binary, defaults to True
    :type useCounter: bool, optional
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
    
    .. _Wildman, Scott A., and Gordon M. Crippen (1999):
        https://pubs.acs.org/doi/abs/10.1021/ci990307l
        
    """
    n_jobs = n_jobs if n_jobs >= 1 else None
    fps = GCfp(useCount)
    pool = Pool(n_jobs)
    fps = pool.map_async(fps.CalculateGCfp, mols).get()
    pool.close()
    pool.join()
    fps = np.array(fps)
    return fps
    
            
def CalculateEState(mols, val=True, n_jobs=1):
    """Atom-based calculation of LogP and MR using Crippen’s approach.
    79 bits
  
    Reference:
        (1) L.B. Kier and L.H. Hall _Molecular Structure Description: The Electrotopological State_ Academic Press (1999)
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param val: If set to True, the fingerprint will presented in the format of vaule of estate, else, would be binary, defauls to True
    :type val: bool, optional
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
    
    """
    n_jobs = n_jobs if n_jobs >= 1 else None
    fps = EStateFP(val)
    pool = Pool(n_jobs)
    fps = pool.map_async(fps.CalculateEState, mols).get()
    pool.close()
    pool.join()
    fps = np.array(fps)
    return fps


def CalculateMACCS(mols, n_jobs=1):
    """There is a SMARTS-based implementation of the 166 public MACCS keys.
    167 bits
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
    
    """
    n_jobs = n_jobs if n_jobs >=1 else None
    fps = MACCS()
    pool = Pool(n_jobs)
    fps = pool.map_async(fps.CalculateMACCS, mols).get()
    pool.close()
    pool.join()
    fps = np.array(fps)
    return fps


def CalculateECFP(mols, radius=2, nBits=1024, n_jobs=1):
    """This family of fingerprints, better known as circular fingerprints, 
    is built by applying the Morgan algorithm to a set of user-supplied atom invariants.
    2^n bits
    
    Reference:
        (1) `Rogers, David, and Mathew Hahn (2010)`_.
    
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param radius: the radius of circle, defaults to 2
    :type rafius: int, optional
    :param nBits: number of bits in the fingerprint, defaults to 1024
    :type nBits: int,optional
    :param useFeatures: to control generate FCFP if True, else ECFP, defaults to False
    :type useFeatures: bool, optional
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
    
    .. _Rogers, David and Mathew Hahn (2010):
        https://pubs.acs.org/doi/abs/10.1021/ci100050t
        
    """
    n_jobs = n_jobs if n_jobs >=1 else None
    fps = Morgan(radius=radius, nBits=nBits)
    pool = Pool(n_jobs)
    fps = pool.map_async(fps.CalculateECFP, mols).get()
    pool.close()
    pool.join()
    fps = np.array(fps)
    return fps


def CalculateFCFP(mols, radius=2, nBits=1024, n_jobs=1):
    """This family of fingerprints, better known as circular fingerprints, 
    is built by applying the Morgan algorithm to a set of user-supplied atom invariants.
    2^n bits
    
    Reference:
        (1) `Rogers, David, and Mathew Hahn (2010)`_.
    
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param radius: the radius of circle, defaults to 2
    :type rafius: int, optional
    :param nBits: number of bits in the fingerprint, defaults to 1024
    :type nBits: int,optional
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
    
    .. _Rogers, David and Mathew Hahn (2010):
        https://pubs.acs.org/doi/abs/10.1021/ci100050t
        
    """
    n_jobs = n_jobs if n_jobs >=1 else None
    fps = Morgan(radius=radius, nBits=nBits)
    pool = Pool(n_jobs)
    fps = pool.map_async(fps.CalculateFCFP, mols).get()
    pool.close()
    pool.join()
    fps = np.array(fps)
    return fps

def CalculateDaylight(mols, minPath=1, maxPath=7, nBits=2048, n_jobs=1):
    """A Daylight-like fingerprint based on hashing molecular subgraphs
    2^n bits
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param minPath: minimum number of bonds to include in the subgraphs, defaults to 1
    :type minPath: int, optional
    :param maxPath: maximum number of bonds to include in the subgraphs, defaults to 7
    :type maxPath: int, optional
    :param nBits: number of bits in the fingerprint, defaults to 2048
    :type nBits: int, optional
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
        
    """
    n_jobs = n_jobs if n_jobs >=1 else None
    fps = Daylight(minPath=minPath, maxPath=maxPath, nBits=nBits)
    pool = Pool(n_jobs)
    fps = pool.map_async(fps.CalculateDaylight, mols).get()
    pool.close()
    pool.join()
    fps = np.array(fps)
    return fps

                                           
def CalculatePubChem(mols, n_jobs=1):
    """Calculate PubChem Fingerprints
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
    
    """
    n_jobs = n_jobs if n_jobs >=1 else None
    fps = PubChem()
    pool = Pool(n_jobs)
    fps = pool.map_async(fps.CalculatePubChem, mols).get()
    pool.close()
    pool.join()
    fps = np.array(fps)
    return fps
 
    
def CalculateIFG(mols, useCount=True, n_jobs=1):
    """An algorithm to identify functional groups in organic molecules
    
    Reference:
        (1) `Peter Ertl (2017)`_.
        
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param useCount: If set to True, the fingerprint will presented in the format of counter(not only 1 and 0) else, would be binary, defaults to True
    :type useCounter: bool, optional
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :return: fingerprints
    :rtype: numpy.ndarray
    
    .. _Peter Ertl (2017):
        https://jcheminf.springeropen.com/articles/10.1186/s13321-017-0225-z
        
    """
    fp = IFG(useCount=useCount, n_jobs=n_jobs)
    fps = fp.CalculateIFG(mols)
    return fps
    
if '__main__' == __name__:
    import time
    import warnings
    warnings.filterwarnings('ignore')
    
    start = time.clock()
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
            ]*100
    
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    fps = CalculateIFG(mols,n_jobs=4)
    end = time.clock()
    print(fps.shape)
    print(end-start)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

