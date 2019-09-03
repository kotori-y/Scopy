# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:53:13 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

♥I love Princess Zelda forever♥
"""

from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem import AllChem as Chem

class Daylight(object):
    """
    Brief:
    -----------
    a Daylight-like fingerprint based on hashing molecular subgraphs
    
    Parameters:
    -----------
    minPath: int(optional, default:1)
        minimum number of bonds to include in the subgraphs
    maxPath: int(optional, default:7)
        maximum number of bonds to include in the subgraphs
    nBits: int(optional, default:2048)
        number of bits in the fingerprint
    """
    def __init__(self, minPath, maxPath, nBits):
        self.minPath = minPath
        self.maxPath = maxPath
        self.nBits = nBits
        
    def CalculateDaylight(self,mol):
        """
        Parameters:
        -----------
        mol: rdkit.Chem.rdchem.Mol
            
        Return:
        -----------
        fp: list
        """
        fp = RDKFingerprint(mol, minPath=self.minPath, maxPath=self.maxPath, fpSize=self.nBits)
        fp = list(fp)
        return fp
    


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
    fp = Daylight(minPath=1, maxPath=7, nBits=2048)
    mol = Chem.MolFromSmiles(smis[3])
    fp = fp.CalculateDaylight(mol)
#    fps = np.array(fps)
    print(fp)

    
    
    
    
    