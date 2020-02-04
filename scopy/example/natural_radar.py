# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:58:28 2020

@author: Kotori Y
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: kotori@cbdd.me& kotori@csu.edu.cn
@Blog: https://blog.iamkotori.com

"""

from multiprocessing import Pool

import pandas as pd
from rdkit import Chem

from scopy.visualize import pc_depict


def GetAndSave(smi):
    """
    Get radar plot consisted of Lipinski and Veber
    
    Rule:
        MW <= 500
        logP <= 5
        nHD <= 5
        nHA <= 10
        TPSA <= 140
    
    :param smi: the SMILES of molecule
    :type smi: string
    :return: plt.figure
    """
    prop_kws = {'MW':[0,500],'logP':[None,5],'nHD':[0,5],
                'nHA':[0,10],'TPSA':[0,140]}
    mol = Chem.MolFromSmiles(smi)
    f = pc_depict.rule_radar(mol,prop_kws)
    # f.savefig('demo.pdf')
    return f



if '__main__' == __name__:
    df = pd.read_csv(r"C:\Users\0720\Desktop\Project\scopy_ref\data\Natural_product\Natural_product-NPS20015.csv")
    smi = df.SMILES[9]
    GetAndSave(smi)