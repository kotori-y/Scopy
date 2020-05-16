# -*- coding: utf-8 -*-

#Created on Mon Sep  2 17:31:14 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me



from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Chem import AllChem as Chem



class Morgan(object):
    """This family of fingerprints, better known as circular fingerprints, 
    is built by applying the Morgan algorithm to a set of user-supplied atom invariants.
    2^n bits
    
    Reference:
        (1) `Rogers, David and Mathew Hahn (2010)`_.
    
    :param radius: the radius of circle, defaults to 2
    :type rafius: int, optional
    :param nBits: number of bits in the fingerprint, defaults to 1024
    :type nBits: int,optional
    :param useFeatures: to control generate FCFP if True, else ECFP, defaults to False
    :type useFeatures: bool, optional
    
    .. _Rogers, David and Mathew Hahn (2010):
        https://pubs.acs.org/doi/abs/10.1021/ci100050t
        
    """
    def __init__(self, radius, nBits):
        """Initialization
        
        """
        self.radius = radius
        self.nBits = nBits
        
    def CalculateECFP(self, mol):
        """Function to compute ECFP fingerprint under useFeatures is True
        
        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol
        :return: fingerprint
        :rtype: list
        
        """
        fp = GetMorganFingerprintAsBitVect(mol, radius=self.radius, nBits=self.nBits)
        fp = list(fp)
        return fp
    
    def CalculateFCFP(self, mol):
        """Function to compute ECFP fingerprint under useFeatures is False
        
        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol
        :return: fingerprint
        :rtype: list
        
        """
        fp = GetMorganFingerprintAsBitVect(mol, radius=self.radius, nBits=self.nBits, useFeatures=True)
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
    fp = Morgan(radius=2,nBits=1024)
    mol = Chem.MolFromSmiles(smis[3])
    fp = fp.CalculateECFP(mol)
#    fps = np.array(fps)
    print(fp)