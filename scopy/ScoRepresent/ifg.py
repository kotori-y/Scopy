# -*- coding: utf-8 -*-

#Created on Thu Dec  5 16:32:11 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me



from multiprocessing import Pool
import warnings, sys
warnings.filterwarnings('ignore')
from collections import Counter
import numpy as np
from rdkit import RDConfig
sys.path.append(RDConfig.RDContribDir)
from IFG.ifg import identify_functional_groups



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


class IFG(object):
    """
    An algorithm to identify functional groups in organic molecules
    Undefined bits(depend to the samples)
    
    :param useCount: If set to True, the fingerprint will presented in the format of counter(not only 1 and 0) else, would be binary, defaults to True
    :type useCounter: bool, optional
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    
    """
    def __init__(self, useCount=False, n_jobs=1):
        """Init
        
        """
        self.useCount = useCount
        self.n_jobs = n_jobs if n_jobs >= 1 else None
        
    def _getifg(self, mols):
        """
        *Internal use only*
        
        :param mols: the molecule to be scanned.
        :type mols: iterable object, each element is rdkit.Chem.rdchem.Mol
        
        """
        pool = Pool(self.n_jobs)
        fgs = pool.map_async(GetIFG, mols).get()
        pool.close()
        pool.join()
        fgs = [[fg['type'] for fg in item] for item in fgs]
        return fgs
    
    def _getfp(self, mols):
        """
        *Internal use only*
        
        :param mols: the molecule to be scanned.
        :type mols: iterable object, each element is rdkit.Chem.rdchem.Mol
        
        """
        fgs = self._getifg(mols)
        total = set(sum(fgs,[]))
        IDX = dict(zip(total, range(len(total))))
        
        for fg in fgs:
            fp = [0]*len(total)
            count = Counter(fg)
            if self.useCount:          
                for k,v in count.items():
                    fp[IDX[k]] = v
            else:
                for k in count.keys():
                    fp[IDX[k]] = 1
            yield np.array(fp)
            
    def CalculateIFG(self, mols):
        """
        Calculate IFG fingrtprint
        
        :param mols: the molecule to be scanned.
        :type mols: iterable object, each element is rdkit.Chem.rdchem.Mol
        :return: the ifg fingerprints
        :rtype: numpy.ndarray
        
        """
        return np.stack(self._getfp(mols))
        
        
        
    
    
if '__main__' == __name__:
    from rdkit import Chem
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
    fp = IFG(n_jobs=4)
    fps = fp.CalculateIFG(mols)
    print(fps.shape)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    