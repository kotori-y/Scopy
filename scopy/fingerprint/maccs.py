# -*- coding: utf-8 -*-

#Created on Thu Nov 28 16:42:49 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me
#
#♥I love Princess Zelda forever♥


from rdkit.Chem.rdMolDescriptors import GetMACCSKeysFingerprint


class MACCS(object):
    """There is a SMARTS-based implementation of the 166 public MACCS keys.
    167 bits
    
    Note:
        (1) Most of the differences have to do with aromaticity
        (2) There’s a discrepancy sometimes because the current RDKit definitions do not require multiple matches to be distinct. e.g. the SMILES C(=O)CC(=O) can match the (hypothetical) key O=CC twice in my definition. It’s not clear to me what the correct behavior is.
        (3) Some keys are not fully defined in the MDL documentation.
        (4) Two keys, 125 and 166, have to be done outside of SMARTS. 
    
    """
    def __init__(self):
        """Initialization
        
        """
        pass
    
    def CalculateMACCS(self, mol):
        """There is a SMARTS-based implementation of the 166 public MACCS keys.
        167 bits
        
        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol
        :return: fingerprint
        :rtype: list
            
        """
        fp = list(GetMACCSKeysFingerprint(mol))
        return fp




if '__main__' == __name__:
    from rdkit.Chem import AllChem as Chem
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

    mol = Chem.MolFromSmiles(smis[0])
    fp = MACCS()
    fp = fp.CalculateMACCS(mol)
    print(fp)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    