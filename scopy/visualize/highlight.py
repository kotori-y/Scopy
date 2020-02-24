# -*- coding: utf-8 -*-

#Created on Mon Sep 16 09:56:26 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me

import re
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import rdDepictor



def HighlightAtoms(mol,highlightAtoms,figsize=[400,200],kekulize=True):
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
    def _revised(svg_words):
        """
        """
        svg_words =  svg_words.replace(
            'stroke-width:2px','stroke-width:1.5px').replace(
                'font-size:17px','font-size:15px').replace(
                    'stroke-linecap:butt','stroke-linecap:square').replace(
                        'fill:#FFFFFF','fill:none').replace(
                        'svg:','')                
        return svg_words
                                                           
    mc = Chem.Mol(mol.ToBinary())
    
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(*figsize)
    drawer.DrawMolecule(mc,highlightAtoms=highlightAtoms)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return SVG(_revised(svg))
    

if '__main__' == __name__:
    mol = Chem.MolFromSmiles('C1=CC=C2C(=O)CC(=O)C2=C1')
    svg = HighlightAtoms(mol, highlightAtoms=(0, 1, 2, 6, 8,10))