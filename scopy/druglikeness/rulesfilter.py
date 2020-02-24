# -*- coding: utf-8 -*-

#Created on Tue Jun 25 21:59:42 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me



from rdkit import Chem
try:
    from . import molproperty
except:
    import sys
    sys.path.append('.')
    import molproperty


def CheckEganRule(mol, detail=False, showSMILES=False):
    """
    Bad or Good oral biovailability rule
    
    Reference:
        Egan, William J., Kenneth M. Merz, and John J. Baldwin. 
        J Med Chem, 43.21 (2000): 3867-3877.
    
    Rule details:
        0 <= TPSA <= 132
        -1 <= logP <=6
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of Egan's rule
    TPSA = molproperty.CalculateTPSA(mol)
    logP = molproperty.CalculateLogP(mol)
    #Determine whether the molecular match each rule
    atPSA = (0 <= TPSA <= 132)
    alogP = (-1 <= logP <= 6)
    #Give the advice
    if atPSA&alogP:
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violsted rules
    nviolate = 2 - (atPSA + alogP)
    #res
    if detail:
        items = ['tPSA', 'logP', 'Disposed','nViolate']
        vals = [TPSA, logP, disposed, nviolate]
    else:
        items = ['Disposed', 'nViolate']
        vals = [disposed, nviolate]
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckVeberRule(mol, detail=False, showSMILES=False):
    """
    Bad or Good oral biovailability rule
    
    Reference:
        Veber, Daniel F., et al.
        Journal of medicinal chemistry 45.12 (2002): 2615-2623.
        
    Rule details:
        nRot <= 10
        TPSA <= 140
        nHB <= 12
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of Veber's rule
    nRot = molproperty.CalculateNumRotatableBonds(mol)
    tPSA = molproperty.CalculateTPSA(mol)
    nHB = molproperty.CalculateNumHyBond(mol)   
    #Determine whether the molecular match each rule
    anRot = (nRot <= 10)
    atPSA = (tPSA <= 140)
    anHB = (nHB <= 12)    
    #Give the advice
    if (anRot&atPSA&anHB):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 3 - (anRot + atPSA + anHB)
    #res
    if detail:
        items = ['nRot', 'tPSA', 'nHB', 'Disposed','nViolate']
        vals = [nRot, tPSA, nHB, disposed, nviolate]
    else:
        items = ['Disposed', 'nViolate']
        vals = [disposed, nviolate]   
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckLipinskiRule(mol, detail=False, showSMILES=False):
    """
    Check molecular under Lipinski's rule
    
    Reference:
        Lipinski, Christopher A., et al.
        Advanced drug delivery reviews 23.1-3 (1997): 3-25.
    
    Rule details:
        MW <= 500(Da)
        logP <= 5
        nHD <= 5
        nHA <= 10
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of Lipinskin's rule
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHD = molproperty.CalculateNumHDonors(mol)
    nHA = molproperty.CalculateNumHAcceptors(mol)   
    #Determine whether the molecular match each rule
    aw = (MW <= 500)
    alogp = (logP <= 5)
    anhd = (nHD <= 5)
    anha = (nHA <= 10)
    #Count the number of matched rules
    nviolate = 4 - (aw + alogp + anhd + anha)
    #Give the disposed
    if nviolate >= 2:
        disposed = 'Rejected'
    else:
        disposed = 'Accepted'
        
    if detail:
        items = ['MW', 'logP', 'nHD', 'nHA', 'Disposed', 'nViolate']
        vals = [MW, logP, nHD, nHA, disposed, nviolate]
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)    
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic   
    

def CheckBeyondRo5(mol, detail=False, showSMILES=False):
    """
    Check molecular under beyond Ro5
    
    Reference:
        Doak, Bradley C., et al.
        journal of medicinal chemistry 59.6 (2015): 2312-2327.
        
    Rule details:
        MW <= 1000
        -2 <= logP <= 10
        nHD <= 6
        nHA <= 15
        PSA <= 250
        nRot <= 20
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of Lipinskin's rule
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHD = molproperty.CalculateNumHDonors(mol)
    nHA = molproperty.CalculateNumHAcceptors(mol)
    PSA = molproperty.CalculateTPSA(mol)
    nRot = molproperty.CalculateNumRotatableBonds(mol)
    #Determine whether the molecular match each rule
    aMW = (MW<=1000)
    alogP = (-2<=logP<=10)
    anHD = (nHD<=6)
    anHA = (nHA<=15)
    aPSA = (PSA<=250)
    aRot = (nRot<=20)
    #Check whether mol pass the whole rules
    if (aMW&alogP&anHD&anHA&aRot&aPSA):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'   
    #Count the number of violated rules
    nviolate = 6 - (aMW + alogP + anHD + anHA + aRot + aPSA)    
    if detail:
        items = ['MW', 'logP', 'nHD', 
                 'nHA', 'PSA','nRot', 
                 'Disposed','nViolate']
        vals = (MW, logP, nHD, nHA, PSA, nRot, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)        
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic   


def CheckPfizerRule(mol, detail=False, showSMILES=False):
    """
    Check molecular under Rfizer Rule(3/75 Rule)
    
    Reference:
        Hughes, Jason D., et al. 
        Bioorganic & medicinal chemistry letters 18.17 (2008): 4872-4875.
        
    Rule details:
        logp > 3
        TPSA < 75
        
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties
    logP = molproperty.CalculateLogP(mol)
    PSA = molproperty.CalculateTPSA(mol)
    #Determine whether the molecular match each rule
    alogP = (logP > 3)
    aPSA = (PSA < 75)    
    #Give the advice
    if alogP&aPSA:
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 2 - (alogP + aPSA)
    #res
    if detail:
        items = ['logP', 'PSA', 
                 'Disposed', 'nViolate']
        vals = (logP,PSA,disposed,nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed,nviolate)    
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic
    

def CheckGSKRule(mol, detail=False, showSMILES=False):
    """
    Check molecular under GSK rule(4/400 Rule)
    
    Reference:
        Gleeson, M. Paul.
        Journal of medicinal chemistry 51.4 (2008): 817-834.
        
    Rule details:
        MW <= 400
        logP <= 4
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)   
    #Determine whether the molecular match each rule    
    aMW = (MW <= 400)
    alogP = (logP <= 4)    
    #Give the advice
    if aMW&alogP:
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 2 - (aMW+alogP)
    #res
    if detail:
        items = ['MW', 'logP', 
                 'Disposed', 'nViolate']
        vals = (MW, logP, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckOralMacrocycles(mol, detail=False, showSMILES=False):
    """
    Check molecular under oral macrocycles rules
    
    Reference:
        Giordanetto, Fabrizio, and Jan Kihlberg.
        Journal of medicinal chemistry 57.2 (2013): 278-295.        

    Rule details:
        MW < 1000
        logP < 10
        nHD < 5
        TPSA < 250
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of Lipinskin's rule
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHD = molproperty.CalculateNumHDonors(mol)
    tPSA = molproperty.CalculateTPSA(mol)   
    #Determine whether the molecular match each rule
    aMW = (MW < 1000)
    alogP = (logP < 10)
    anHD = (nHD < 5)
    aPSA = (tPSA < 250)  
    #Give the advice
    if (aMW&alogP&anHD&aPSA):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violate rules
    nviolate = 4 - (aMW + alogP + anHD + aPSA)
    #res
    if detail:
        items = ['MW', 'logP', 'nHD', 'tPSA', 
                 'Disposed', 'nViolate']
        vals = (MW, logP, nHD, tPSA, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate) 
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckOpreaRule(mol, detail=False, showSMILES=False):
    """
    Reference:
        Oprea, Tudor I.
        Journal of computer-aided molecular design 14.3 (2000): 251-264.
        
    Rules details:
        nRing >= 3
        nRig >= 18
        nRot >=6
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of Oprea
    nRing = molproperty.CalculateNumRing(mol)
    nRig = molproperty.CalculateNumRigidBonds(mol)
    nRot = molproperty.CalculateNumRotatableBonds(mol)
    #Determine whether the molecular match each rule
    anRing = (nRing >= 3)
    anRig = (nRig >= 18)
    anRot = (nRot >= 6)
    #Give the advice
    if (anRing&anRig&anRot):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 3 - (anRing + anRig + anRot)
    #Res
    if detail:
        items = ['nRing', 'nRig', 'nRot', 
                 'Disposed','nViolate']
        vals = (nRing, nRig, nRot, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic

    

# def CheckOpreaTwoRule(mol):
#     """
#     Reference.:
               
#     Rule details:
#         mw <= 450
#         -3.5 <= clogP <= 4.5
#         -4 <= logD <= 4
#         nring <= 4
#         no. of nonterminal single bonds<=10 #this property should be realized
#         NHD <= 5
#         NHA <= 8
#     """
#     pass 


def CheckGhoseRule(mol, detail=False, showSMILES=False):
    """
    Check molecular under Ghose rule
    
    Reference.:
        Ghose, Arup K., Vellarkad N. Viswanadhan, and John J. Wendoloski. 
        Journal of combinatorial chemistry 1.1 (1999): 55-68.
    
    Rules details:
        -0.4 < logP < 5.6
        160 < MW < 480
        40 < MR< 130
        20 < nAtom < 70
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of Oprea
    logP = molproperty.CalculateLogP(mol)
    MW = molproperty.CalculateMolWeight(mol)
    MR = molproperty.CalculateMolMR(mol)
    nAtom = molproperty.CalculateNumAtoms(mol)    
    #Determine whether the molecular match each rule
    alogP = (-0.4 < logP < 5.6)
    aMW = (160 < MW < 480)
    aMR = (40 < MR< 130)
    anAtom = (20 < nAtom<  70)
    #Give the advice
    if (alogP&aMW&aMR&anAtom):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 4 - (alogP + aMW + aMR + anAtom)
    #Res
    if detail:
        items = ['logP', 'MW', 'MR', 'nAtom', 
                 'Disposed', 'nViolate']
        vals = (logP, MW, MR, nAtom, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


# def CheckKelderRule(mol, detail=False, showSMILES=False):
#     """
#     Check moleculars under Kelder rules
    
#     Ref.:
#         Kelder, Jan, et al.
#         Pharm Res, 16.10 (1999): 1514-1519.
        
#     Rule details:
#         tpsa <120 for orally active;
#         tpsa < 60-70 for brain penetration
        
#     """
#     pass


def CheckREOS(mol, detail=False, showSMILES=False):
    """
    Check molecular under REOS program
    
    Reference:
        Walters, W. Patrick, and Mark Namchuk.
        Nat Rev Drug Discov, 2.4 (2003): 259.
        
    Rule details:
        200 <= MW <= 500
        -5 <= logP <= 5
        nHD <= 5
        nHA <= 10
        nRot <= 8
        TPSA <= 150
        -4 <= fChar <= 4
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of REOS
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHD = molproperty.CalculateNumHDonors(mol)
    nHA = molproperty.CalculateNumHAcceptors(mol)
    nRot = molproperty.CalculateNumRotatableBonds(mol)
    TPSA = molproperty.CalculateTPSA(mol)
    fChar = molproperty.CalculateMolFCharge(mol) 
    #Determine whether the molecular match each rule
    amw = (200 <= MW <= 500)
    alogP = (-5 <= logP <= 5)
    anhd = (nHD <= 5)
    anha = (nHA <= 10)
    anrot = (nRot <= 8)
    atpsa = (TPSA <= 150)
    afchar = (-4 <= fChar <= 4)
    #Give the advice
    if (amw&alogP&anhd&anha&anrot&atpsa&afchar):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 7 - (amw + alogP + anhd + anha + anrot + atpsa + afchar)
    #Res
    if detail:
        items = ['MW', 'logP', 'nHD', 'nHA',
                 'nRot', 'TPSA', 'fChar', 
                 'Disposed','nViolate']
        vals = (MW, logP, nHD, nHA, nRot, TPSA, fChar, disposed,nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckGoldenTriangle(mol, detail=False, showSMILES=False):
    """
    Check molecular under 'Golden Triangle'
    
    Reference:
        Johnson, Ted W., Klaus R. Dress, and Martin Edwards.
        Bioorg Med Chem Lett, 19.19 (2009): 5560-5564.
        
    Rule details:
        200 <= MW <= 500
        -2 <= logD <= 5
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of REOS
    MW = molproperty.CalculateMolWeight(mol)
    logD = molproperty.CalculateLogD(mol)
    #Determine whether the molecular match each rule
    amw = (200 <= MW <= 500)
    alogd = (-2 <= logD <=5)
    #Give the advice
    if (amw&alogd):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 2 - (amw + alogd)
    #Res
    if detail:
        items = ['MW', 'logD', 
                 'Disposed', 'nViolate']
        vals = (MW, logD, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckXuRule(mol, detail=False, showSMILES=False):
    """
    Check molecular under Xu's rule
    
    Reference:
        
    
    Rule details:
        nhd <= 5
        nha <= 10
        3 <= rot <= 35
        1 <= nring <= 7
        10 <= nhev <= 50
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of Xu's rule
    nHD = molproperty.CalculateNumHDonors(mol)
    nHA = molproperty.CalculateNumHAcceptors(mol)
    nRot = molproperty.CalculateNumRotatableBonds(mol)
    nRing = molproperty.CalculateNumRing(mol)
    nHev = molproperty.CalculateNumHeavyAtom(mol)   
    #Determine whether the molecular match each rule
    anHD = (nHD <= 5)
    anHA = (nHA <= 10)
    anRot = (3 <= nRot <= 35)
    anRing = (1 <= nRing <= 7)
    anHev = (10 <= nHev <= 50)
    #Give the advice
    if (anHD&anHA&anRot&anRing&anHev):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 5 - (anHD + anHA + anRot + anRing + anHev)
    #res
    if detail:
        items = ['nHD', 'nHA', 'nRot', 
                 'nRing', 'nHev', 'Disposed', 'nViolate']
        vals = (nHD, nHA, nRot, nRing, nHev, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


# def CheckSchneiderRule(mol):
#     """
#     Check  molecular under Schneider rule
    
#     Reference:
#         Schneider, Nadine, et al.
#         J Chem Inf Model, 48.3 (2008): 613-628.
        
#     Rule details:
#         mw > 230
#         nhd > 0
#         nha > 0
#         nrot > 0
#         nring > 0
#         mr > 40
#         functional groups > 0
#         molecular volume > 191
#     """
#     pass


def CheckRo4(mol, detail=False, showSMILES=False):
    """
    Referenece:
        
    Rule details:
        MW <= 400
        logP <= 4
        nHD <= 4
        nHA <= 8
        TPSA <= 120
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of REOS
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHD = molproperty.CalculateNumHDonors(mol)
    nHA = molproperty.CalculateNumHAcceptors(mol)
    tPSA = molproperty.CalculateTPSA(mol)    
    #Determine whether the molecular match each rule
    aMW = (MW <= 400)
    alogP = (logP <= 4)
    anHD = (nHD <= 4)
    anHA = (nHA <= 8)
    atPSA = (tPSA <= 120)   
    #Give the advice
    if (aMW&alogP&anHD&anHA&atPSA):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 5 - (anHD + anHA + aMW + alogP + atPSA)
    #res
    if detail:
        items = ['MW', 'logP', 'nHD', 
                 'nHA', 'tPSA', 'Disposed', 'nViolate']
        vals = (MW, logP, nHD, nHA, tPSA ,disposed, nviolate)
    else:
        items = ['Disposed','nViolate']
        vals = (disposed,nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckRo3(mol, detail=False, showSMILES=False):
    """
    Check molecular under Ro3
    
    Ref.:
        Congreve, Miles, et al.
        Drug discovery today 19.8 (2003): 876-877.
        
    Rule details:
        MW <= 300
        -3 <= logP <= 3
        NHD <= 3
        NHA <= 6
        PSA <= 60
        nRot <= 3
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of REOS
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHD = molproperty.CalculateNumHDonors(mol)
    nHA = molproperty.CalculateNumHAcceptors(mol)
    tPSA = molproperty.CalculateTPSA(mol)
    nRot = molproperty.CalculateNumRotatableBonds(mol)    
    #Determine whether the molecular match each rule
    aMW = (MW <= 300)
    alogP = (-3 <= logP <= 3)
    anHD = (nHD <= 3)
    anHA = (nHA <= 6)
    atPSA = (tPSA <= 60)
    aRot = (nRot <= 3)
    #Give the advice
    if (aMW&alogP&anHD&anHA&atPSA&anRot):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 6 - (aMW + alogP + anHD + anHA + atPSA + anRot)
    #res
    if detail:
        items = ['MW', 'logP', 'nHD', 
                 'nHA','tPSA','nRot', 
                 'Disposed','nViolate']
        vals = (MW, logP, nHD, nHA, tPSA, nRot, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)    
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckRo2(mol, detail=False, showSMILES=False):
    """
    Check molecular under RO2
    
    Ref.:
        Goldberg, Frederick W., et al.
        Drug Discovery Today 20.1 (2015): 11-17.
        
    Rule details:
        MW <= 200
        Logp <= 2
        NHD <= 2
        NHA <= 4

    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of REOS
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHD = molproperty.CalculateNumHDonors(mol)
    nHA = molproperty.CalculateNumHAcceptors(mol)    
    #Determine whether the molecular match each rule
    aMW = (MW <= 200)
    alogP = (logP <= 2)
    anHD = (nHD <= 2)
    anHA = (nHA <= 4)    
    #Give the advice
    if (aMW&alogP&anHD&anHA):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 4 - (aMW + alogP + anHD + anHA)
    #res
    if detail:
        items = ['MW', 'logP', 'nHD', 
                 'nHA', 'Disposed', 'nViolate']
        vals = (MW, logP, nHD, nHA, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


# def DrugLikeOne(mol):
#     """
#     Based on Drug-Like Soft designed by FAFDrugs4(http://fafdrugs4.mti.univ-paris-diderot.fr/index.html)
    
#     Through combining several articles describing drugs' physico-chemical properties 
#     and an in-house statistical analysis of drugs.
    
#     Reference:
#         (1) Lipinski, Christopher A., et al.
#         Adv Drug Deliv Rev, 223.1-3 (1997): 3-25.
#         (2) Oprea, Tudor I.
#         J Comput Aided Mol Des, 14.3 (2000): 251-264.
#         (3) Irwin, John J., and Brian K. Shoichet.
#         J Chem Inf Model, 45.1 (2005): 177-182.
#         (4) Oprea, Tudor I., et al.
#         J Chem Inf Comput Sci, 41.5 (2001): 1308-1315.
#         (5) Pihan, Emilie, et al.
#         Bioinformatics, 28.11 (2012): 1540-1541.
        
#     Rule details:
#         100 <= weight <= 600
#         -3 <= logP <= 6
#         nha <= 12
#         nhd <= 7
#         TPSA <= 180
#         nrot <= 11
#         nrig <= 30
#         nring <= 6
#         maxring <= 18
#         3 <= ncarb <=35
#         1 <= nhet <= 15
#         0.1 <= hetcar <=1.1
#         ncharged <= 4
#         -4 <= totalchar <= 4
#     """
#     pass


# def DrugLikeTwo(mol):
#     """
#     Based on Drug-Like Soft designed by FAFDrugs4(http://fafdrugs4.mti.univ-paris-diderot.fr/index.html)
    
#     Through combining several articles describing drugs' physico-chemical properties 
#     and an in-house statistical analysis of drugs.
    
#     Reference.:
#         (1) Oprea, Tudor I., et al.
#         J Comput Aided Mol Des, 14.3 (2000): 251-264.
#         (2) Workman, Paul, and Ian Collins.
#         Chem Biol, 17.6 (2010): 561-577.
#         (3) Baell, Jonathan B.
#         J Chem Inf Model, 53.1 (2012): 39-55.
#         (4) Brenk, Ruth, et al.
#         ChemMedChem 3.3 (2008): 435-444.
        
#     Rule details:
#         150 <= weight <= 400
#         -3 <= logP <= 4
#         nha <= 7
#         nhd <= 4
#         TPSA <= 160
#         nrot <= 9
#         nrig <= 30
#         nring <= 4
#         maxring <= 18
#         3 <= ncarb <=35
#         1 <= nhet <= 15
#         0.1 <= hetcar <=1.1
#         ncharged <= 4
#         -4 <= totalchar <= 4
#         nstero <= 2
#     """
#     pass


# def CheckZinc(mol):
#     """
#     Check mol under ZINC
    
#     Reference.:
#         Irwin, John J., and Brian K. Shoichet. 
#         J Chem Inf Model, 45.1 (2005): 177-182.
        
#     Rule details:
#         60 <= weight <= 100
#         -4 <= logP <= 6
#         nha <= 11
#         nhd <= 6
#         TPSA <= 150
#         nrot <= 12
#         nrig <= 50
#         nring <= 7
#         maxring <= 12
#         ncarb >= 3
#         nhet >= 0
#         hetcar <= 2.0
#         ncharged <= 4
#         -4 <= totalchar <= 4
#         nstero <= 2
#     """
#     pass


def CheckCNS(mol, detail=False, showSMILES=False):
    """
    Check mol under CNS
    
    Reference:
        Jeffrey, Phil, and Scott Summerfield.
        Neurobiol Dis, 37.1 (2010): 33-37.
        
    Rule details:
        135 <= weight <= 582
        -0.2 <= logP <= 6.1
        nha <= 5
        nhd <= 3
        3 <= TPSA <= 118
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of REOS
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHD = molproperty.CalculateNumHDonors(mol)
    nHA = molproperty.CalculateNumHAcceptors(mol)
    tPSA = molproperty.CalculateTPSA(mol)
    #Determine whether the molecular match each rule
    aMW = (135 <= MW <= 582)
    alogP = (-0.2 <= logP <= 6.1)
    anHD = (nHD <= 3)
    anHA = (nHA <= 5)  
    atPSA = (3 <= tPSA <= 118)
    #Give the advice
    if (aMW&alogP&anHD&anHA&atPSA):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 5 - (aMW + alogP + anHD + anHA + atPSA)
    #res
    if detail:
        items = ['MW', 'logP', 'nHD', 
                 'nHA', 'tPSA', 'Disposed', 'nViolate']
        vals = (MW, logP, nHD, nHA, tPSA, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def CheckRespiratory(mol, detail=False, showSMILES=False):
    """
    Check mol under Respiratory
    
    Reference:
        Ritchie, Timothy J., Christopher N. Luscombe, and Simon JF Macdonald. 
        J Chem Inf Model, 49.4 (2009): 1025-1032.
        
    Rule details:
        240<=MW<= 520
        -2.0<=logP<=4.7
        6<=nHB<=12
        51<=tPSA<=135
        3<=nRot<=8
        1<=nRing<=5
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: bool, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    #Calculate required properties of REOS
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)
    nHB = molproperty.CalculateNumHyBond(mol)
    tPSA = molproperty.CalculateTPSA(mol)
    nRot = molproperty.CalculateNumRotatableBonds(mol)
    nRing = molproperty.CalculateNumRing(mol)
    #Determine whether the molecular match each rule
    aMW = (240 <= MW <= 520)
    alogP = (-0.2 <= logP <= 4.7)
    anHB = (6 <= nHB <= 12)  
    atPSA = (51 <= tPSA <= 135)
    anRot = (3 <= nRot <= 8)
    anRing = (1 <= nRing <= 5)
    #Give the advice
    if (aMW&alogP&anHB&atPSA&anRot&anRing):
        disposed = 'Accepted'
    else:
        disposed = 'Rejected'
    #Count the number of violated rules
    nviolate = 6 - (aMW + alogP + anHB + atPSA + anRot + anRing)
    #res
    if detail:
        items = ['MW', 'logP', 'nHB', 
                 'tPSA', 'nRot', 'nRing', 
                 'Disposed', 'nViolate']
        vals = (MW, logP, nHB, tPSA, nRot, nRing, disposed, nviolate)
    else:
        items = ['Disposed', 'nViolate']
        vals = (disposed, nviolate)
    
    dic = {'SMILES':Chem.MolToSmiles(mol)} if showSMILES else {}
    dic.update(dict(zip(items, vals)))
    return dic


def Check_CustomizeRule(mol, prop_kws, closed_interval=True, detail=False, showSMILES=False):
    """
    You could customize the rule with mostly properties ypu want.
    
    :param mol: molecule
    :type mol: rdkit.Chem.rdchem.Mol
    :param prop_kws: the keys of dict are properties you want to check; the values should be a tuple or list with two elements,present the left- and right-bounded respectively.
    :type prop_kws: `dict`
    :param closed_interval: control whether using closed interval, defaults to True
    :type closed_interval: `bool`, optional
    :param detail: Control returning specific infomation or not, defaults to False
    :type detail: `bool`, optional
    :param showSMILES: Control returning SMILES or not, defaults to False
    :type showSMILES: bool, optional
    
    :return: Result after scanning. If detail has been set to False, only return 'Disposed' and 'nViolate', otherwise exrta returning each property 
    :rtype: `dict`
    
    """
    allowed = ['MW','Vol','Dense','fChar','nBond','nAtom','nHD','nHA','nHB',
             'nHet','nStero','nHev','nRot','nRig','nRing',
             'logP','logD','pKa','logSw','ab','MR','tPSA','AP','HetRatio',
             'Fsp3','MaxRing','QEDmean','QEDmax','QEDnone','SAscore','NPscore',
             'nSingle','nDouble','nTriple','nC','nB','nF','nCl','nBr','nI',
             'nP','nS','nO','nN']
    
    properties = molproperty.GetProperties(mol)
    keys = list(prop_kws.keys())
    try:
        res = [properties[key] for key in keys]
        if closed_interval:
            bo = [((not prop_kws.get(key)[0] or properties[key] >= prop_kws.get(key)[0])and\
                   (not prop_kws.get(key)[-1] or properties[key] <= prop_kws.get(key)[-1]))\
                   for key in keys]
        else:
            bo = [((not prop_kws.get(key)[0] or properties[key] > prop_kws.get(key)[0])and\
                   (not prop_kws.get(key)[-1] or properties[key] < prop_kws.get(key)[-1]))\
                   for key in keys]
    except KeyError:
        raise ValueError('the key of prop_kws must a member of {}'.format(', '.join(allowed)))
    
    VioProp = [keys[idx] for idx,item in enumerate(bo) if not item]
    nViolate = len(keys)-sum(bo)
    
    if detail:       
        keys.extend(['nViolate','VioProp'])
        res.extend([nViolate,VioProp]) 
    else:
        keys = ['nViolate','VioProp']
        res = [nViolate,VioProp]
    
    dic = dict(zip(keys, res))
    if showSMILES:
        dic['SMILES'] = Chem.MolToSmiles(mol)
    return dic
       


if __name__ =='__main__':
    
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]','CC(=O)OC1=CC=CC=C1C(=O)O']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
    smiring = ['C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2','C1CCC2C3CC(C4CCCC5CCCCC45)CCC3CCC2C1','C1CCC2(CCC3(CCCCC3)CC2)CC1']
    for index, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        print('Index:{}'.format(index))
#        res = Check_CustomizeRule(mol, prop_kws={'MW':[0,100],
#                                                 'nN':[0,2]},detail=False, showSMILES=False)
        res = CheckGoldenTriangle(mol)
        print(res)
    
    
    
    
    
    
    
    
    
    
    
    
    
   