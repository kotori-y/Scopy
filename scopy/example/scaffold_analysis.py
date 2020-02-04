# -*- coding: utf-8 -*-

#Created on Wed Jan  1 17:35:09 2020
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me
#
#♥I love Princess Zelda forever♥


import csv
from rdkit import Chem
import openbabel as ob
import pandas as pd #This package should be installed
from scopy.structure_alert.SmartsFilter import Filter
from scopy.druglikeness.druglikeness import PC_properties, PC_rules
from scopy.visualize import mcloud
import warnings
warnings.filterwarnings('ignore')
import datetime
import os
from multiprocessing import Pool



def obsmitosmile(smi):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "can")
    conv.SetOptions("K", conv.OUTOPTIONS)
    mol = ob.OBMol()
    conv.ReadString(mol, smi)
    smile = conv.WriteString(mol)
    smile = smile.replace('\t\n', '')
    return smile


def getmol(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        return mol
    else:
        return Chem.MolFromSmiles(obsmitosmile(smi))


def main(folder,file):
    """
    """
    os.chdir(folder)
    
    print('+++++++++++++PREPARE++++++++++++++++')
    
    data = pd.read_csv(file)
#    data['mol'] = data.mol.map(lambda x: obsmitosmile(x))
    name,_ = os.path.splitext(file)
        
    ps = Pool(24)
    mols = ps.map_async(getmol, data.mol.values).get()
    ps.close()
    ps.join()
    
    scount_muc = mcloud.CountScaffold(mols)
    scount_cbr = mcloud.CountScaffold(mols,stype='Carbon')
    
    with open(r'C:\Users\0720\Desktop\station\scopy\doc3\image\scaffolds.txt', 'w', newline='') as f_obj:
        writer = csv.writer(f_obj, delimiter='\t')
        for k,v in scount_muc.items():
             writer.writerow([k,v])
    f_obj.close()
    
    with open(r'C:\Users\0720\Desktop\station\scopy\doc3\image\scaffolds.txt', 'w', newline='') as f_obj:
        writer = csv.writer(f_obj, delimiter='\t')
        for k,v in scount_cbr.items():
             writer.writerow([k,v])
    f_obj.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    