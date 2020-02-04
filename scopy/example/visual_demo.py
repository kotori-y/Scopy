# -*- coding: utf-8 -*-

#Created on Sat Jan 11 15:01:43 2020
#
#@Author: Zhi-Jiang Yang, Jie Dong, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.iamkotori.com
#
#♥I love Princess Zelda forever♥

import os

import pandas as pd

import openbabel as ob
from rdkit import Chem
from scopy.visualize import highlight


def obsmitosmile(smi):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "can")
    conv.SetOptions("K", conv.OUTOPTIONS)
    mol = ob.OBMol()
    conv.ReadString(mol, smi)
    smile = conv.WriteString(mol)
    smile = smile.replace('\t\n', '')
    return smile


def getMol(smi):
    m = Chem.MolFromSmiles(smi)
    if m:
        return m
    else:
        m = Chem.MolFromSmiles(obsmitosmile(smi))
        return m


def ShowSubstructure(file):
    """
    
    """
    res = pd.read_csv(file)
    res = res[res.Disposed=='Rejected']
    
    res['MatchedAtoms'] = res.MatchedAtoms.map(lambda x: eval(x))
    mols = (getMol(smi) for smi in res.SMILES.values)
    m_atoms,m_names = res.MatchedAtoms.values,res.MatchedNames.values
    
    for mol,atoms,names in zip(mols,m_atoms,m_names):
        for atom in atoms:
            idx = 0
            for highl in atom:
                idx += 1
                fig = highlight.HighlightAtoms(mol,highl)
                return fig
            
            
            
if '__main__' == __name__:
    fig = ShowSubstructure(r"C:\Users\0720\Desktop\Project\scopy_ref\data\Withdraw_drug\Res\Genotoxic_Carcinogenicity_Mutagenicity.csv")
    
    
    
    
    
    
    
    
    
    
    
    
    
    