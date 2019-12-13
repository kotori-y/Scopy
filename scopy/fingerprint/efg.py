# -*- coding: utf-8 -*-

#Created on Tue Sep  3 10:31:34 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me
#
#♥I love Princess Zelda forever♥



import pickle, gzip, os, csv
try:
    from ..structure_alert.SmartProcess import _CheckPattl
    from .. import ScoConfig
except:
    import sys
    sys.path.append('..')
    from structure_alert.SmartProcess import _CheckPattl
    import ScoConfig
    
from rdkit import Chem


def _Generatepkl(endpoint):
    """
    *Internal Use Only*
    
    the pkl file in this package, storing the rdkit.Chem.rdchem.Mol object,
    was generated under the environment whose rdkit version is '2018.03.4'.
    Since, the file may can't be successfully loaded. This function is designed for
    re-generating a pkl file under this situation.
    
    :param endpoint: the name of file
    :return: None
    
    """
    file = os.path.join(ScoConfig.EFGDir,'{}'.txt.format(endpoint))
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
    out = pickle.dumps(lines,protocol=-1)
    outfile = os.path.join(ScoConfig.EFGDir,'{}.pkl.gz'.format(endpoint))
    with gzip.open(outfile,'wb') as f_out:
        f_out.write(out)
    f_out.close()


class EFG(object):
    """classification system termed “extended functional groups” (EFG),
    which are an extension of a set previously used by the CheckMol software, 
    that covers in addition heterocyclic compound classes and periodic table groups. 
    583 bits
    
    Reference:
        (1) `Salmina, Elena, Norbert Haider and Igor Tetko (2016)`_.
    
    :param useCount: If set to True, the fingerprint will presented in the format of counter(not only 1 and 0) else, would be binary, defaults to True
    :type useCounter: bool, optional
    
    .. _Salmina, Elena, Norbert Haider and Igor Tetko (2016):
        https://www.mdpi.com/1420-3049/21/1/1
    
    """
    def __init__(self, useCount):
        """Initialization
        
        """
        self.useCount = useCount
        file = os.path.join(ScoConfig.EFGDir,'Extended_Functional_Groups.pkl.gz')
        try:
            pattl = pickle.load(gzip.open(file,'rb'))
        except:
            _Generatepkl('Extended_Functional_Groups')
            pattl = pickle.load(gzip.open(file,'rb'))    
        self.rejected = [patt[-2] for patt in pattl]
        self.accepted = [patt[-1] for patt in pattl]      
        self.name = [patt[0] for patt in pattl]
        
    
    def CalculateEFG(self,mol):
        """Calulating EFG fingerprint
        
        Note:
            following bits are always to be set as 1 or 0: bit0, bit1, bit2, bit59, bit60, bit61, bit62, bit209, bit218, bit544, bit545.
            It's diffcult to counter them like others bits, because these bits contain "not"
        
        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol
        :return: fingerprint
        :rtype: list
        
        """
        Bool = _CheckPattl(mol, self.rejected, self.accepted)
        if not self.useCount:
            fp = [x+0 for x in Bool]
        else:
            fp = [0]*583
            for idx,bo,rejl,accl in zip(range(583),Bool,self.rejected,self.accepted):
                if bo:
                    if not rejl:
                        fp[idx] = sum([len(mol.GetSubstructMatches(acc)) for acc in accl])
                    else:
                        fp[idx] = 1
                else:
                    pass
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
    fp = EFG(useCount=0)
    mol = Chem.MolFromSmiles(smis[3])
    fp = fp.CalculateEFG(mol)
#    fps = np.array(fps)
    print(fp)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    