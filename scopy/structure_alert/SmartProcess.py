# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 21:59:42 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

import _pickle as cPickle
import gzip
from rdkit.Chem import AllChem as Chem
from collections import namedtuple
try:
	from .. import ScoConfig
except:
	import sys
	sys.path.append('..')
	import ScoConfig

__doc__ = """
Some SMARTSs in Ochem database(https://ochem.eu/alerts/home.do) would 
contain logic symbol, include "AND", "OR" and "NOT", cannot be recognized 
by RDKit. So that we split all SMARTS to two parts, "Rejected" and "Accepted",
meant should not be matched and be necessarily matched part respectively 
according to the SMARTS.

e.g.
-----------
>>> SMARTS: [CX4]F AND NOT [$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]
    - > Reject: [$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]
    - > Accept: [CX4]F
    It means if a molecule don't match "[$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]"
    and match [CX4]F will be considered as matching SMARTS: "[CX4]F AND NOT [$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]"

>>> SMARTS: NOT [!#1!#6]
    Reject: [!#1!#6]
    Accept: None
    It means if only a molecule don't match "[!#1!#6]" will be considered as matching "NOT [!#1!#6]"

>>> SMARTS: [SX2][Sv6X4](=[OX1])(=[OX1])
    Reject: None
    Accept: [SX2][Sv6X4](=[OX1])(=[OX1])
    It means if only a molecule match "[SX2][Sv6X4](=[OX1])(=[OX1])" will be considered as matching "[SX2][Sv6X4](=[OX1])(=[OX1])"
"""


def _Generatepkl(endpoint):
    """  
    *Internal Use Only*
    
    the pkl file in this package, storing the rdkit.Chem.rdchem.Mol object,
    was generated under the environment whose rdkit version is '2018.03.4'.
    Since, the file may can't be successfully loaded. This function is designed for
    re-generating a pkl file under this situation.
    
    Parameters:
    -----------
    endpoint: str
        the name of file
        
    Return:
    -----------
    None
    """
    import csv 
    import os
    file = os.path.join(ScoConfig.SmartDir,'{}'.txt.format(endpoint))
    with open(file,'r',encoding='utf-8') as f_obj:
        lines = csv.reader(f_obj,delimiter='\t')
        next(lines)  
        lines = tuple(lines)
    f_obj.close()
    
    for line in lines:
        rej,acc = line[-2],line[-1]
        if rej:
            rej = eval(rej) 
            rej = [Chem.MolFromSmarts(x) for x in rej]
            line[-2] = rej         
        if acc:
            acc = eval(acc) 
            acc = [Chem.MolFromSmarts(x) for x in acc]
            line[-1] = acc          
    out = cPickle.dumps(lines,protocol=-1)
    outfile = os.path.join(ScoConfig.PattDir,'{}.pkl.gz'.format(endpoint))
    with gzip.open(outfile,'wb') as f_out:
        f_out.write(out)
    f_out.close()
    
        
def _Loadpkl(endpoint):
    """ 
    *Internal Use Only*
    
    loading the specific pkl file which contain the 'Rejected' and 'Accepted' SMARTS
    in rdkit.Chem.rdchem.Mol object format to avoid repeated reading SMARTS by 'MolFromSmarts'
    
    Parameters:
    -----------
    endpoint: string
        the endpoint of pkl file meant
    
    Return:
    -----------
    pattl: list
        whose element ia also a list with four elements:
        0: the name of SMARTS, 1:original SMARTS,
        2: the 'rejected' part of SMARTS,
        3: the 'accepted' part of SMARTS.      
    """
    filename = ScoConfig.PattDir +'\\{}.pkl.gz'.format(endpoint)
    try:
        pattl = cPickle.load(gzip.open(filename,'rb'))
    except FileNotFoundError:
        class EndpointError(Exception):
            pass
        import os
        files = os.listdir(ScoConfig.SmartDir)
        raise EndpointError('Endpoint must be one of these:\n\t{}.'\
                            .format(', '.join([file.replace('.txt','') for file in files])))
    except:
        _Generatepkl(endpoint)
        _Loadpkl(endpoint)
    return pattl


def _CheckPattl(mol, rejected_pattl, accepted_pattl):
    """
    *Internal Use Only*
    
    Checking mol through 'rejected' and 'accepted' part respectively
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
        the molecule to be scanned
    reject_pattl: rdkit.Chem.rdchem.Mol
        the 'rejected' part of SMARTS
    accepted_pattl: rdkit.Chem.rdchem.Mol
        the 'accepted' part of SMARTS
        
    Yield:
    -----------
    bool, True meant matched SMARTS, else unmatched.
    """    
    for reject_list,accept_list in zip(rejected_pattl,accepted_pattl):
        if (not reject_list) or (not any([mol.HasSubstructMatch(patt) for patt in reject_list])):
            if accept_list:
                yield all([mol.HasSubstructMatch(patt) for patt in accept_list])
            else:
                yield True
        else:
            yield False
      

def _CheckWithSmarts(mol, pattl, endpoint, detail=False):
    """  
    *Internal Use Only*
    
    checking molecule(s) wheather or not 
    has(have) some (toxic) substructure(s) through comparing the 
    SMILES of molecule(s) and the SMARTS of (toxic) substructure(s) 
    which obtained from Ochem Database(https://ochem.eu/alerts/home.do).
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
        the molecule to be scanned
    pattl: list
        generated by _Loadpkl, containing 'rejected' and 'accepted' part
    endpoint: string
        the endpoint name of pattl that de used to scan mol
    detail: bool(default: False)
        if True, the function would return more specific infomation
        Warning: If true, would be time and space consumed.
    
    Return:
    -----------
    A namedtuple:
        if detail, would return 'Disposed','MatchedAtoms','MatchedNames','Endpoint':
            Disposed: bool
                true meant unmatched any SMARTS in this endpoint, else at least one to be matched
            MatchedAtoms: list
                which atoms has matched SMARTS
            MatchedNames: list
                which SMARTS to be matched
            Endpoint: string
                the name of endpoint 
        else only return 'Disposed' and 'Endpoint'
    """ 
    reject_pattl = map(lambda x: x[-2], pattl)
    accept_pattl = map(lambda x: x[-1], pattl)
      
    if detail:
        reject_pattl = tuple(reject_pattl)
        accept_pattl = tuple(accept_pattl)
                
        matched_names = []
        matched_atoms = []
        
        names = map(lambda x: x[0], pattl)

        temp = _CheckPattl(mol,reject_pattl,accept_pattl)
        temp = tuple(temp)
        if any(temp):               
            disposed = False
            matched_patts = (pattl for pattl,bo in zip(accept_pattl,temp) if bo)
            for patts in matched_patts:
                if patts:
                    matched_sub = [mol.GetSubstructMatches(patt) for patt in patts]
                else:
                    matched_sub = [(tuple([atom.GetIdx() for atom in mol.GetAtoms()]),)]
                matched_atoms.extend(matched_sub) 
            matched_names = list((name for name,bo in zip(names,temp) if bo))       
       
        else:
            disposed = True
            matched_atoms = ['-']
            matched_names = ['-']
                        
  
        res = namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint'])
        del reject_pattl,accept_pattl
        return res(Disposed=disposed,MatchedAtoms=matched_atoms,MatchedNames=matched_names,Endpoint=endpoint)
    
    else:
        temp = _CheckPattl(mol,reject_pattl,accept_pattl)
        if not any(temp):
            disposed = True
        else: 
            disposed = False
        res = namedtuple('CheckRes',['Disposed','Endpoint'])
        return res(Disposed=disposed,Endpoint=endpoint)  
                

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
    mol = Chem.MolFromSmiles(smis[0])
    pattl = _Loadpkl(endpoint='Pains')
    res = _CheckWithSmarts(mol,pattl=pattl,endpoint='Pains',detail=1)
    print(res)
#    pic = VisualizeFragment(mol,atoms=[0,2,1])
