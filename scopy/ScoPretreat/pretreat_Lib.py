# -*- coding: utf-8 -*-
"""
Created on Sat May 23 23:09:14 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from multiprocessing import Pool
from rdkit import Chem
from .pretreat import StandardMol



def StandardMols(mols, n_jobs=1):
    """
    """
    mols = mols if type(mols) is not Chem.rdchem.Mol else [mols]
    n_jobs = n_jobs if n_jobs>=1 else None
    
    pool = Pool(n_jobs)
    sm = pool.map_async(StandardMol, mols).get()
    pool.close()
    pool.join()
    
    return sm




if '__main__' == __name__:
    smis = ['O=C([O-])c1ccccc1','C[n+]1c([N-](C))cccc1','[2H]C(Cl)(Cl)Cl']
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    sm = StandardMols(mols, n_jobs=3)
    stsmis = [Chem.MolToSmiles(mol) for mol in sm]
    print(stsmis)
    