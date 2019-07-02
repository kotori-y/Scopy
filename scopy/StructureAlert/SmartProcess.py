# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 21:59:42 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

import csv
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from collections import namedtuple
from scopy import ScoConfig

__doc__ = """

Some SMARTS in Ochem contain contain logic symbol, include "AND", "OR" and "NOT" 
cannot be recognized by RDKit. So that we split all SMARTS to two parts, "Reject" and "Accept",
meant shouldn't be marched and necessarily be matched part respectively.

---
e.g.
>>> SMARTS: [CX4]F AND NOT [$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]
Reject: [$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]
Accept: [CX4]F
It means if a molecule don't match "[$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]"
and match [CX4]F will be considered as matching SMARTS: "[CX4]F AND NOT [$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]"

>>> NOT [!#1!#6]
Reject: [!#1!#6]
Accept: None
It means if only a molecule don't match "[!#1!#6]" will be considered as matching "NOT [!#1!#6]"

>>> [SX2][Sv6X4](=[OX1])(=[OX1])
Reject: None
Accept: [SX2][Sv6X4](=[OX1])(=[OX1])
It means if only a molecule match "[SX2][Sv6X4](=[OX1])(=[OX1])" will be considered as matching "[SX2][Sv6X4](=[OX1])(=[OX1])"

"""



def _GetPattl(smartlist):
    """
    """
    if smartlist == '':
        return None    
    else:
        pattl = [Chem.MolFromSmarts(x) for x in eval(smartlist)]
        return pattl


def _CheckPattl(mol,reject_pattl,accept_pattl):
    """
    """
    res = []
    for reject_list,accept_list in zip(reject_pattl,accept_pattl):
        if reject_list is None or (reject_list is not None\
                                   and not any([mol.HasSubstructMatch(patt) for patt in reject_list])):
            if accept_list is not None:
                res.append(all([mol.HasSubstructMatch(patt) for patt in accept_list]))
            else:
                res.append(True)
        else:
            res.append(False)
    return res


def _CheckWithSmarts(mol,endpoint='all',show=False):
    """
    #################################################################
    *Internal Use Only*
    
    checking molecule(s) wheather or not 
    has(have) some (toxic) substructure(s) through comparing the 
    SMILES of molecule(s) and the SMARTS of (toxic) substructure(s) 
    which obtained from ToxAlerts Database(https://ochem.eu/alerts/home.do).
    #################################################################
    """
    filename = ScoConfig.SmartDir +'\\{}.txt'.format(endpoint)
    matched_names = []
    matched_atoms = []
    
    with open(filename,'r',encoding='utf-8') as f_obj:
        
        lines = csv.reader(f_obj,delimiter='\t')
        next(lines)  
        lines = tuple(lines)
        
        names = map(lambda x: x[0], lines)
        names = tuple(names)
        
        reject_smarts = tuple(map(lambda x: x[-2], lines))
        reject_pattl = tuple(map(_GetPattl,reject_smarts))
       
        accept_smarts = tuple(map(lambda x: x[-1], lines))
        accept_pattl = tuple(map(_GetPattl,accept_smarts))

        temp = _CheckPattl(mol,reject_pattl,accept_pattl)
        if endpoint !='Extended_Functional_Groups':  
            if any(temp):
                disposed = 'Rejected'
                matched_patts = [accept_pattl[index]\
                                 for index,bo in enumerate(temp) if bo is True]
                for patts in matched_patts:
                    if patts is not None:
                        matched_sub = [mol.GetSubstructMatches(patt) for patt in patts]
                    else:
                        matched_sub = [(tuple([atom.GetIdx() for atom in mol.GetAtoms()]),)]
                    matched_atoms.extend(matched_sub)
                matched_names = [names[index]\
                                 for index,bo in enumerate(temp) if bo is True]        
           
            else:
                disposed = 'Accepted'
                matched_atoms = ['-']
                matched_names = ['-']
                                
            if show:
                res = namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint'])
                checkres = res(Disposed=disposed,MatchedAtoms=matched_atoms,MatchedNames=matched_names,Endpoint=endpoint)
            else:
                res = namedtuple('CheckRes',['Disposed','Endpoint'])
                checkres = res(Disposed=disposed,Endpoint=endpoint)
        else:
            checkres = ''.join([['0','1'][x] for x in temp])
    
    f_obj.close()    
    return checkres


def VisualizeFragment(mol,atoms):
    """
    #################################################################
    This function is used for show which part of fragment matched the SMARTS
    
    -Usage:
        
        pic = Visualize(mol,atoms)
        
        Input: mol is a molecule object. atoms is a tuple
        
        Output pic is a PIL.Image.Image
    #################################################################
    """
    
    pic = Draw.MolToImage(mol,highlightAtoms=atoms,size=(600,600),fitImage=True)
    return pic




if '__main__' == __name__:
    mol = Chem.MolFromSmiles('C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2') 
    res = _CheckWithSmarts(mol,endpoint='Extended_Functional_Groups')
#    pic = VisualizeFragment(mol,atoms=[0,2,1])
    print(res)