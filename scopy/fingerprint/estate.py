# -*- coding: utf-8 -*-

#Created on Wed Jul 17 10:15:43 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me



from rdkit.Chem.EState import EState
from rdkit.Chem.EState import Fingerprinter as ESFP
from rdkit import Chem

Version=1.0
# This module is derived from our previous work

################################################################

class EStateFP(object):
    """
    79 bits
    
    Reference:
        (1) L.B. Kier and L.H. Hall _Molecular Structure Description: The Electrotopological State Academic Press (1999)
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    :param val: If set to True, the fingerprint will presented in the format of vaule of estate, else, would be binary, defauls to True
    :type val: bool, optional
    
    """
    def __init__(self, val=True):
        """Initialization
        
        """
        self.val = val
    
    def CalculateEState(self, mol):
        """Calculate EState fingerprint
        
        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol
        :return: fingerprint
        :rtype: list
        
        """
        if self.val:
            fp = ESFP.FingerprintMol(mol)[1]
        else:
            fp = ESFP.FingerprintMol(mol)[0]
        return fp

def _CalculateEState(mol, skipH=1):
    """
    **Internal used only**
    
    Get the EState value of each atom in a molecule
    """
    res = EState.EStateIndices(mol)
    return res


if __name__=='__main__':
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
    fp = EStateFP(val=True)
    mol = Chem.MolFromSmiles(smis[3])
    fp = fp.CalculateEState(mol)
#    fps = np.array(fps)
    print(fp)


    
    
    
    
    
    
    
    
      
    
    
    