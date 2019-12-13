# -*- coding: utf-8 -*-

#Created on Mon Sep 16 09:56:26 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me


from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor
    

def HighlightAtoms(mol,highlightAtoms,figsize=[400,200]):
    """This function is used for showing which part of fragment matched the SMARTS by the id of atoms.
    
    :param mol: The molecule to be visualized
    :type mol: rdkit.Chem.rdchem.Mol
    :param highlightAtoms: The atoms to be highlighted
    :type highlightAtoms: tuple
    :param figsize: The resolution ratio of figure
    :type figsize: list
    :return: a figure with highlighted molecule
    :rtype: IPython.core.display.SVG
    
    """    
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(*figsize)
    drawer.DrawMolecule(mol,highlightAtoms=highlightAtoms)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    fig = SVG(svg)
    return fig


if '__main__' == __name__:
    from rdkit import Chem
    mol = Chem.MolFromSmiles('C1=CC=C2C(=O)CC(=O)C2=C1')
    svg = HighlightAtoms(mol, highlightAtoms=(0, 1, 2, 6, 8,10))