# -*- coding: utf-8 -*-

#Created on Mon Sep  2 09:24:41 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me



import sys
try:
    from .. import ScoConfig
except:
    sys.path.append('..')
    import ScoConfig
import csv
from rdkit.Chem import AllChem as Chem

class GCfp(object):
    """Atom-based calculation of LogP and MR using Crippen’s approach.
    110 bits
    
    Reference:
        (1) `Wildman, Scott A., and Gordon M. Crippen (1999)`_.
    
    :param useCount: If set to True, the fingerprint will presented in the format of counter(not only 1 and 0) else, would be binary, defaults to True
    :type useCounter: bool, optional
    
    .. _Wildman, Scott A., and Gordon M. Crippen (1999):
        https://pubs.acs.org/doi/abs/10.1021/ci990307l
        
    """   
    def __init__(self, useCount=True):
        """Initialization
        
        """
        with open(ScoConfig.CrippenDir + '\\Crippen.txt') as f_obj:
            lines = csv.reader(f_obj,delimiter='\t')
            next(lines)
            self.pattl = tuple(map(lambda line: Chem.MolFromSmarts(line[1]), lines))
        f_obj.close()        
        self.useCount = useCount

    def CalculateGCfp(self,mol):
        """Calculate GC fingerprint
        
        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol
        :return: fingerprint
        :rtype: list
        
        """
        nAtom = len(mol.GetAtoms())
        doneAtoms = [0]*nAtom
        fp = [0]*110
        idx = -1
        flag = 0  
        found = 0
        
        for patt in self.pattl:
            idx += 1
            count = 0
            for match in mol.GetSubstructMatches(patt):
                firstidx = match[0]
                if not doneAtoms[firstidx]:
                    doneAtoms[firstidx] = 1
                    count += 1
                    found += 1
                    if found >= nAtom:
                        flag = 1
                        break
            fp[idx] = count
            if flag:
                if not self.useCount:
                    fp = [0 if not int(bit) else 1 for bit in fp]
                return fp
        
        else:
            return [0]*110


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

    fp = GCfp(useCount=False)
    mol = Chem.MolFromSmiles(smis[3])
    fp = fp.CalculateGCfp(mol)
#    fps = np.array(fps)
    print(fp)























