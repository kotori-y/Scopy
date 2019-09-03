# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 13:43:48 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

from scopy.StructureAlert import SmartsFilter
from rdkit import Chem
import pandas as pd
import timeit
import multiprocessing

df = pd.read_csv(r"C:\Users\0720\Desktop\MATE\yzy\chemblock_r.csv")
df = df.sample(n=10000,random_state=10)
smis = list(df.smiles)

def singleCore_single():
    for smi in smis:
        mol = Chem.MolFromSmiles(smi)
        Hts = SmartsFilter.HTS(mol)
        Hts.scan()
 
def singleCore_library():
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    print('over')
    Hts = SmartsFilter.HTS(mol=mols, stype='library')
    Hts.scan()
#    res = Hts.ScanResult
#    return res
    

def cut(a,b):
    smi = smis[a:b]
    mols = [Chem.MolFromSmiles(s) for s in smi]
    Hts = SmartsFilter.HTS(mol=mols, stype='library')
    Hts.scan()
    
        
        
if '__main__' == __name__:   
    print(len(smis))     
#    print('>>>SingleCore under stype=single timecost: {}'.format(timeit.timeit(stmt=singleCore_single,number=3)),end='\n\n')
    print('>>>SingleCore under stype=library timecost: {}'.format(timeit.repeat(stmt=singleCore_library,number=3,repeat=2)),end='\n\n')
    
    def multiCore():
        pool = multiprocessing.Pool()
        pool.apply_async(cut,args=(0,251))
        pool.apply_async(cut,args=(251,501))
        pool.apply_async(cut,args=(501,751))
        pool.apply_async(cut,args=(751,1001))
        pool.close()
        pool.join()
    
    print('>>>MultiCore timecost: {}'.format(timeit.repeat(stmt=multiCore,number=3,repeat=2)))

#
#     