# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:55:57 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from functools import partial
from multiprocessing import Pool
from collections import Counter
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold



def _getscaffold(mol,stype='Murcko'):
    """
    *Internal used only*
    
    """
    assert stype in ['Murcko','Carbon'
                     ], 'scaffold type must be a member of "Murcko" or "Carbon"'
    core = MurckoScaffold.GetScaffoldForMol(mol)
    core = core if stype=='Murcko' else MurckoScaffold.MakeScaffoldGeneric(core)
    return Chem.MolToSmiles(core, isomericSmiles=False, canonical=True)


def CountMurckoFramework(mols):
    """
    Counting the frequency of Murcko Framework
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    
    
    :return: the SMILES of framework and its frequency
    :rtype: dict
    
    """
    fn = partial(_getscaffold,stype='Murcko')
    ps = Pool(4)
    scaffolds = ps.map_async(fn, mols).get()
    ps.close()
    ps.join()
    count = dict(Counter(scaffolds))
    try:
        del count['']
    except KeyError:
        pass
    
    return count


def CountCarbonScaffold(mols):
    """
    Counting the frequency of Carbon Scaffold
    
    :param mols: the molecule to be scanned.
    :type mols: Iterable object, each element is rdkit.Chem.rdchem.Mol
    
    
    :return: the SMILES of framework and its frequency
    :rtype: dict
    
    """
    fn = partial(_getscaffold,stype='Carbon')
    ps = Pool(4)
    scaffolds = ps.map_async(fn, mols).get()
    ps.close()
    ps.join()
    count = dict(Counter(scaffolds))
    try:
        del count['']
    except KeyError:
        pass
    
    return count



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
            'CCCCC(=O)[H]', #Biodegradable
            'C1=CN=C(C(=O)O)C=C1', #Chelating
            'C(OC)1=CC=C2OCC3OC4C=C(OC)C=CC=4C(=O)C3C2=C1',
            'C1=C2N=CC=NC2=C2N=CNC2=C1', #Genotoxic_Carcinogenicity_Mutagenicity
            'N(CC)(CCCCC)C(=S)N', #Idiosyncratic
            ]
    
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    scount = CountCarbonScaffold(mols)
    print(scount)






