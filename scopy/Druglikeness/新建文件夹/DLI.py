# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 12:36:00 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

from scopy.Druglikeness import CalculateProperty


def rule_1(mol):
    """
    number of non-H atoms
    """
    return CalculateProperty.CalculateNumHeavyAtom(mol)


def rule_2(mol):
    """
    The number of SSSR
    """
    return CalculateProperty.CalculateNumRing(mol)


def rule_3(mol):
    """
    the molecular cyclized degree (MCD)
    Eq.: MCD = 100 - 100*(nHev/(nBond+))
    """
    nHev = rule_1(mol)
    nBond = CalculateProperty.CalculateNumBonds(mol)    
    nRing = CalculateProperty.CalculateNumRing(mol)
    if nRing == 0:
        return 0
    else:
        ringinfo = mol.GetRingInfo()
        rings = ringinfo.AtomRings()
        cores = [list(rings[0])]
        for ring in rings[1:]:
            core = []
        pass


    