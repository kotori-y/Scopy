# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 21:59:42 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

__doc__ = """
    This module implments properties obtained from module molproperty
    ---
    We have collected following rules:

    Egan Rule     0<=tPSA<=132; -1<=logP<=6
    Veber Rule   nRot<= 10; tPSA<= 140; nHB<= 12
    LipinskiRule  MW<=500; logP<=5, nHD<=5, nHA <=10
    BeyondRo5   Mw<=1000; -2<=logP<=10; nHD<=6, nHA<=15; tPSA<=250; nRot<=20
    Pfizer Rule     logP>3; PSA<75
    GSK Rule        MW<=400; logP<=4
    OralMacrocycles     MW<1000; logP<10; nHD<5; PSA<250
    Oprea Rule   nRing>=3,nRig>=18,nRot>=6  
    Ghose Rule      -0.4<logP<5.6; 160<MW<480; 40<MR<130; 20<nAtom<70
    Xu Rule     nHD<=5; nHA <= 10; 3 <= rot <= 35; 1 <= nring <= 7; 10 <= nhev <= 50
    Ro4 Rule    MW<=400; logP<=4; nHD<=4; NHA<=8; PSA<=120
    Ro3 Rule    MW<=300; -3<=logP<=3; nHD<=3; nHA<=6; PSA<=60
    Ro2 Rule    MW<=200; logP<=2; nHD<=2; nHA<=4
    REOS Rule    200<=MW<=500; -5<=logP<=5; nHD<=5; nHA<=10; nRot<=8; TPSA<=150; -4<=fChar<=4
    GoldenTriangle
    ---
    Followed should be achieved in the future:

    OpreaTwo Rule: mw <= 450; -3.5<=clogP<=4.5; -4<=logD<=4; nring<= 4; no. of nonterminal single bonds<=10; nHD<=5; nHA<=8
    Kelder Rule
    Schneider Rule
    DrugLikeOne
    DrugLikeTwo
    Zinc
    Cns
    Respiratory
    """

try:
    from . import molproperty
except:
    import sys
    sys.path.append('.')
    import molproperty
from collections import namedtuple
from rdkit import Chem


def CheckEganRule(mol, detail=False):
    """
    Bad or Good oral biovailability rule
    
    Ref.:
    -----------
    Egan, William J., Kenneth M. Merz, and John J. Baldwin. 
    J Med Chem, 43.21 (2000): 3867-3877.
    
    Rule details:
    -----------
    0 <= TPSA <= 132
    -1 <= logP <=6
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            bRo5Ruler(MW=180.16, logP=1.31, nHD=1, nHA=3, PSA=63.6, nRot=3, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
    """
    #Calculate required properties of Lipinskin's rule
    TPSA = molproperty.CalculateTPSA(mol)
    logP = molproperty.CalculateLogP(mol)
    #Determine whether the molecular match each rule
    atPSA = (0 <= TPSA <= 132)
    alogP = (-1 <= logP <= 6)
    #Give the advice
    if atPSA&alogP:
        disposed = True
    else:
        disposed = False
    #Count the number of violsted rules
    violate = 2 - (atPSA + alogP)
    #res
    if detail:
        res = namedtuple('EganRuler', ['tPSA', 'logP', 
                                       'Disposed','nViolate'])
        checkres = res(TPSA, logP, disposed, violate)
    else:
        res = namedtuple('EganRuler', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def CheckVeberRule(mol, detail=False):
    """
    Bad or Good oral biovailability rule
    
    Ref.:
    -----------
    Veber, Daniel F., et al.
    Journal of medicinal chemistry 45.12 (2002): 2615-2623.
        
    Rule details:
    -----------
    nRot <= 10
    TPSA <= 140
    nHB <= 12
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            VeberRuler(nRot=3, tPSA=63.6, nHB=4, Disposed=True, Violate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
    """
    #Calculate required properties of Lipinskin's rule
    nRot = molproperty.CalculateNumRotatableBonds(mol)
    tPSA = molproperty.CalculateTPSA(mol)
    nHB = molproperty.CalculateNumHyBond(mol)   
    #Determine whether the molecular match each rule
    anRot = (nRot <= 10)
    atPSA = (tPSA <= 140)
    anHB = (nHB <= 12)    
    #Give the advice
    if (anRot&atPSA&anHB):
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 3 - (anRot + atPSA + anHB)
    #res
    if detail:
        res = namedtuple('VeberRuler', ['nRot', 'tPSA', 'nHB', 
                                        'Disposed','nViolate'])
        checkres = res(nRot, tPSA, nHB, disposed, violate)
    else:
        res = namedtuple('VeberRuler', ['Disposed','nViolate'])
        checkres = res(disposed, violate)    
    return checkres


def CheckLipinskiRule(mol, detail=False):
    """
    Check molecular under Lipinski's rule
    
    Ref.:
    -----------
    Lipinski, Christopher A., et al.
    Advanced drug delivery reviews 23.1-3 (1997): 3-25.
    
    Rule details:
    -----------
    MW <= 500(Da)
    logP <= 5
    nHD <= 5
    nHA <= 10
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            LipinskiRule(MW=180.16, logP=1.31, nHD=1, nHA=3, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = False
    else:
        disposed = True
        
    if detail:
        res = namedtuple('LipinskiRule', ['MW', 'logP', 'nHD', 
                                          'nHA', 'Disposed', 'nViolate'])
        checkres = res(MW, logP, nHD, nHA, disposed, nviolate)
    else:
        res = namedtuple('LipinskiRule', ['Disposed', 'nViolate'])
        checkres = res(disposed, nviolate)    
    return checkres   
    

def CheckBeyondRo5(mol, detail=False):
    """
    Check molecular under beyond Ro5
    
    Ref.:
    -----------
    Doak, Bradley C., et al.
    journal of medicinal chemistry 59.6 (2015): 2312-2327.
        
    Rule details:
    -----------
    MW <= 1000
    -2 <= logP <= 10
    nHD <= 6
    nHA <= 15
    PSA <= 250
    nRot <= 20
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            bRo5Ruler(MW=180.16, logP=1.31, nHD=1, nHA=3, PSA=63.6, nRot=3, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = True
    else:
        disposed = False   
    #Count the number of violated rules
    violate = 6 - (aMW + alogP + anHD + anHA + aRot + aPSA)    
    if detail:
        res = namedtuple('bRo5Ruler', ['MW', 'logP', 'nHD', 
                                       'nHA', 'PSA','nRot', 
                                       'Disposed','nViolate'])
        checkres = res(MW, logP, nHD, nHA, PSA, nRot, disposed, violate)
    else:
        res = namedtuple('bRo5Ruler', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)        
    return checkres    


def CheckPfizerRule(mol, detail=False):
    """
    Check molecular under Rfizer Rule(3/75 Rule)
    
    Ref.:
    -----------
    Hughes, Jason D., et al. 
    Bioorganic & medicinal chemistry letters 18.17 (2008): 4872-4875.
        
    Rule details:
    -----------
    logp > 3
    TPSA < 75
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            PfizerRule(logP=1.31, PSA=63.6, Disposed=False, nViolate=1)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
    """
    #Calculate required properties
    logP = molproperty.CalculateLogP(mol)
    PSA = molproperty.CalculateTPSA(mol)
    
    #Determine whether the molecular match each rule
    alogP = (logP > 3)
    aPSA = (PSA < 75)    
    #Give the advice
    if alogP&aPSA:
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 2 - (alogP + aPSA)
    #res
    if detail:
        res = namedtuple('PfizerRule', ['logP', 'PSA', 
                                        'Disposed', 'nViolate'])
        checkres = res(logP,PSA,disposed,violate)
    else:
        res = namedtuple('PfizerRule', ['Disposed', 'nViolate'])
        checkres = res(disposed,violate)    
    return checkres
    

def CheckGSKRule(mol, detail=False):
    """
    Check molecular under GSK rule(4/400 Rule)
    
    Ref.:
    -----------
    Gleeson, M. Paul.
    Journal of medicinal chemistry 51.4 (2008): 817-834.
        
    Rule details:
    -----------
    MW <= 400
    logP <= 4
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            GSKRuler(MW=180.16, logP=1.31, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
    """
    #Calculate required properties
    MW = molproperty.CalculateMolWeight(mol)
    logP = molproperty.CalculateLogP(mol)   
    #Determine whether the molecular match each rule    
    aMW = (MW <= 400)
    alogP = (logP <= 4)    
    #Give the advice
    if aMW&alogP:
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 2 - (aMW+alogP)
    #res
    if detail:
        res = namedtuple('GSKRuler', ['MW', 'logP', 
                                      'Disposed', 'nViolate'])
        checkres = res(MW, logP, disposed, violate)
    else:
        res = namedtuple('GSKRuler', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def CheckOralMacrocycles(mol, detail=False):
    """
    Check molecular under oral macrocycles rules
    
    Ref.:
    -----------
    Giordanetto, Fabrizio, and Jan Kihlberg.
    Journal of medicinal chemistry 57.2 (2013): 278-295.        

    Rule details:
    -----------
    MW < 1000(Da)
    logP < 10
    nHD < 5
    TPSA < 250
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            OralMacrocycles(MW=180.16, logP=1.31, nHD=1, tPSA=63.6, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = True
    else:
        disposed = False
    #Count the number of violate rules
    violate = 4 - (aMW + alogP + anHD + aPSA)
    #res
    if detail:
        res = namedtuple('OralMacrocycles', ['MW', 'logP', 'nHD', 'tPSA', 
                                             'Disposed', 'nViolate'])
        checkres = res(MW, logP, nHD, tPSA, disposed, violate)
    else:
        res = namedtuple('OralMacrocycles', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate) 
    return checkres


def CheckOpreaRule(mol, detail=False):
    """
    Ref.:
    -----------
    Oprea, Tudor I.
    Journal of computer-aided molecular design 14.3 (2000): 251-264.
        
    Rules details:
    -----------
    nRing >= 3
    nRig >= 18
    nRot >=6
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            OpreaRule(nRing=1, nRig=8, nRot=3, Disposed=False, nViolate=3)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 3 - (anRing + anRig + anRot)
    #Res
    if detail:
        res = namedtuple('OpreaRule', ['nRing', 'nRig', 'nRot', 
                                       'Disposed','nViolate'])
        checkres = res(nRing, nRig, nRot, disposed, violate)
    else:
        res = namedtuple('OpreaRule', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres

    

def CheckOpreaTwoRule(mol):
    """
    #################################################################
    -Ref.:
               
    -Rule details:
        mw <= 450
        -3.5 <= clogP <= 4.5
        -4 <= logD <= 4
        nring <= 4
        no. of nonterminal single bonds<=10 #this property should be realized
        NHD <= 5
        NHA <= 8
    #################################################################
    """
    pass 


def CheckGhoseRule(mol, detail=False):
    """
    Check molecular under Ghose rule
    
    Ref.:
    -----------
    Ghose, Arup K., Vellarkad N. Viswanadhan, and John J. Wendoloski. 
    Journal of combinatorial chemistry 1.1 (1999): 55-68.
    
    Rules details:
    -----------
    -0.4 < logP < 5.6
    160 < MW < 480
    40 < MR< 130
    20 < nAtom < 70
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            GhoseRule(logP=1.31, MW=180.16, MR=44.71, nAtom=21, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 4 - (alogP + aMW + aMR + anAtom)
    #Res
    if detail:
        res = namedtuple('GhoseRule', ['logP', 'MW', 'MR', 'nAtom', 
                                      'Disposed', 'nViolate'])
        checkres = res(logP, MW, MR, nAtom, disposed, violate)
    else:
        res = namedtuple('GhoseRule', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def CheckKelderRule(mol, detail=False):
    """
    #################################################################
    Check moleculars under Kelder rules
    
    -Ref.:
        Kelder, Jan, et al.
        Pharm Res, 16.10 (1999): 1514-1519.
        
    -Rule details:
        tpsa <120 for orally active;
        tpsa < 60-70 for brain penetration
        
    #################################################################
    """
    pass
    
    return None


def CheckREOS(mol, detail=False):
    """
    Check molecular under REOS program
    
    Ref.:
    -----------
    Walters, W. Patrick, and Mark Namchuk.
    Nat Rev Drug Discov, 2.4 (2003): 259.
        
    Rule details:
    -----------
    200 <= MW <= 500
    -5 <= logP <= 5
    nHD <= 5
    nHA <= 10
    nRot <= 8
    TPSA <= 150
    -4 <= fChar <= 4
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            REOSRule(MW=180.16, logP=1.31, nHD=1, nHA=3, nRot=3, TPSA=63.6, fChar=0, Disposed=False, nViolate=1)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 7 - (amw + alogP + anhd + anha + anrot + atpsa + afchar)
    #Res
    if detail:
        res = namedtuple('REOSRule', ['MW', 'logP', 'nHD', 'nHA',
                                     'nRot', 'TPSA', 'fChar', 
                                     'Disposed','nViolate'])
        checkres = res(MW, logP, nHD, nHA, nRot, TPSA, fChar, disposed,violate)
    else:
        res = namedtuple('REOSRule', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def CheckGoldenTriangle(mol, detail=False):
    """
    Check molecular under 'Golden Triangle'
    
    Ref.:
    -----------
    Johnson, Ted W., Klaus R. Dress, and Martin Edwards.
    Bioorg Med Chem Lett, 19.19 (2009): 5560-5564.
        
    Rule details:
    -----------
    200 <= MW <= 500
    -2 <= logD <= 5
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            GoldenTriangle(MW=180.16, logD=0.6117536115454374, Disposed=False, nViolate=1)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
    """
    #Calculate required properties of REOS
    MW = molproperty.CalculateMolWeight(mol)
    logD = molproperty.CalculateLogD(mol)
    #Determine whether the molecular match each rule
    amw = (200 <= MW <= 500)
    alogd = (-2 <= logD <=5)
    #Give the advice
    if (amw&alogd):
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 2 - (amw + alogd)
    #Res
    if detail:
        res = namedtuple('GoldenTriangle', ['MW', 'logD', 
                                            'Disposed', 'nViolate'])
        checkres = res(MW, logD, disposed, violate)
    else:
        res = namedtuple('GoldenTriangle', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def CheckXuRule(mol, detail=False):
    """
    Check molecular under Xu's rule
    
    Ref.:
    -----------
        
    
    Rule details:
    -----------
    nhd <= 5
    nha <= 10
    3 <= rot <= 35
    1 <= nring <= 7
    10 <= nhev <= 50
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            XuRule(nHD=1, nHA=3, nRot=3, nRing=1, nHev=13, Disposed=True, Violate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 5 - (anHD + anHA + anRot + anRing + anHev)
    #res
    if detail:
        res = namedtuple('XuRule', ['nHD', 'nHA', 'nRot', 
                                    'nRing', 'nHev', 'Disposed', 'nViolate'])
        checkres = res(nHD, nHA, nRot, nRing, nHev, disposed, violate)
    else:
        res = namedtuple('XuRule', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def CheckSchneiderRule(mol):
    """
    Check  molecular under Schneider rule
    
    Ref.:
    -----------
    Schneider, Nadine, et al.
    J Chem Inf Model, 48.3 (2008): 613-628.
        
    Rule details:
        mw > 230
        nhd > 0
        nha > 0
        nrot > 0
        nring > 0
        mr > 40
        ######################
        Followed should be realized
        
        functional groups > 0
        molecular volume > 191
        
        ######################
    #################################################################
    """
    pass
    return None


def CheckRo4(mol, detail=False):
    """
    Ref.:
    -----------
    
    
    Rule details:
    -----------
    MW <= 400
    logP <= 4
    nHD <= 4
    nHA <= 8
    TPSA <= 120
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            Ro4(MW=180.16, logP=1.31, nHD=1, nHA=3, tPSA=63.6, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively. 
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
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 5 - (anHD + anHA + aMW + alogP + atPSA)
    #res
    if detail:
        res = namedtuple('Ro4', ['MW', 'logP', 'nHD', 
                                 'nHA', 'tPSA', 'Disposed', 'nViolate'])
        checkres = res(MW, logP, nHD, nHA, tPSA ,disposed, violate)
    else:
        res = namedtuple('Ro4',['Disposed','nViolate'])
        checkres = res(disposed,violate)
    return checkres


def CheckRo3(mol, detail=False):
    """
    Check molecular under Ro3
    
    Ref.:
    -----------
    Congreve, Miles, et al.
    Drug discovery today 19.8 (2003): 876-877.
        
    Rule details:
    -----------
    MW <= 300
    -3 <= logP <= 3
    NHD <= 3
    NHA <= 6
    PSA <= 60
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            Ro3(MW=180.16, logP=1.31, nHD=1, nHA=3, tPSA=63.6, nRot=3, Disposed=False, nViolate=1)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively. 
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
    anRot = (nRot <= 3)
    #Give the advice
    if (aMW&alogP&anHD&anHA&atPSA&anRot):
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 6 - (aMW + alogP + anHD + anHA + atPSA + anRot)
    #res
    if detail:
        res = namedtuple('Ro3', ['MW', 'logP', 'nHD', 
                                 'nHA','tPSA','nRot', 
                                 'Disposed','nViolate'])
        checkres = res(MW, logP, nHD, nHA, tPSA, nRot, disposed, violate)
    else:
        res = namedtuple('Ro3', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)    
    return checkres


def CheckRo2(mol, detail=False):
    """
    Check molecular under RO2
    
    Ref.:
    -----------
        Goldberg, Frederick W., et al.
        Drug Discovery Today 20.1 (2015): 11-17.
        
    Rule details:
    -----------
    MW <= 200
    Logp <= 2
    NHD <= 2
    NHA <= 4

    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            Ro2(MW=180.16, logP=1.31, nHD=1, nHA=3, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.     
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
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 4 - (aMW + alogP + anHD + anHA)
    #res
    if detail:
        res = namedtuple('Ro2', ['MW', 'logP', 'nHD', 
                                 'nHA', 'Disposed', 'nViolate'])
        checkres = res(MW, logP, nHD, nHA, disposed, violate)
    else:
        res = namedtuple('Ro2', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def DrugLikeOne(mol):
    """
    #################################################################
    Based on Drug-Like Soft designed by FAFDrugs4(http://fafdrugs4.mti.univ-paris-diderot.fr/index.html)
    
    Through combining several articles describing drugs' physico-chemical properties 
    and an in-house statistical analysis of drugs.
    
    -Ref.:
        (1) Lipinski, Christopher A., et al.
        Adv Drug Deliv Rev, 223.1-3 (1997): 3-25.
        (2) Oprea, Tudor I.
        J Comput Aided Mol Des, 14.3 (2000): 251-264.
        (3) Irwin, John J., and Brian K. Shoichet.
        J Chem Inf Model, 45.1 (2005): 177-182.
        (4) Oprea, Tudor I., et al.
        J Chem Inf Comput Sci, 41.5 (2001): 1308-1315.
        (5) Pihan, Emilie, et al.
        Bioinformatics, 28.11 (2012): 1540-1541.
        
    -Rule details:
        100 <= weight <= 600
        -3 <= logP <= 6
        nha <= 12
        nhd <= 7
        TPSA <= 180
        nrot <= 11
        nrig <= 30
        nring <= 6
        maxring <= 18
        3 <= ncarb <=35
        1 <= nhet <= 15
        0.1 <= hetcar <=1.1
        ncharged <= 4
        -4 <= totalchar <= 4
    #################################################################
    """
    pass



def DrugLikeTwo(mol):
    """
    #################################################################
    Based on Drug-Like Soft designed by FAFDrugs4(http://fafdrugs4.mti.univ-paris-diderot.fr/index.html)
    
    Through combining several articles describing drugs' physico-chemical properties 
    and an in-house statistical analysis of drugs.
    
    -Ref.:
        (1) Oprea, Tudor I., et al.
        J Comput Aided Mol Des, 14.3 (2000): 251-264.
        (2) Workman, Paul, and Ian Collins.
        Chem Biol, 17.6 (2010): 561-577.
        (3) Baell, Jonathan B.
        J Chem Inf Model, 53.1 (2012): 39-55.
        (4) Brenk, Ruth, et al.
        ChemMedChem 3.3 (2008): 435-444.
        
    -Rule details:
        150 <= weight <= 400
        -3 <= logP <= 4
        nha <= 7
        nhd <= 4
        TPSA <= 160
        nrot <= 9
        nrig <= 30
        nring <= 4
        maxring <= 18
        3 <= ncarb <=35
        1 <= nhet <= 15
        0.1 <= hetcar <=1.1
        ncharged <= 4
        -4 <= totalchar <= 4
        nstero <= 2
    #################################################################
    """
    pass


def CheckZinc(mol):
    """
    #################################################################
    Check mol under ZINC
    
    -Ref.:
        Irwin, John J., and Brian K. Shoichet. 
        J Chem Inf Model, 45.1 (2005): 177-182.
        
    -Rule details:
        60 <= weight <= 100
        -4 <= logP <= 6
        nha <= 11
        nhd <= 6
        TPSA <= 150
        nrot <= 12
        nrig <= 50
        nring <= 7
        maxring <= 12
        ncarb >= 3
        nhet >= 0
        hetcar <= 2.0
        ncharged <= 4
        -4 <= totalchar <= 4
        nstero <= 2
    #################################################################
    """
    pass
    
    return None


def CheckCNS(mol, detail=False):
    """
    Check mol under CNS
    
    Ref.:
    -----------
    Jeffrey, Phil, and Scott Summerfield.
    Neurobiol Dis, 37.1 (2010): 33-37.
        
    Rule details:
    -----------
    135 <= weight <= 582
    -0.2 <= logP <= 6.1
    nha <= 5
    nhd <= 3
    3 <= TPSA <= 118
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            CNS(MW=180.16, logP=1.31, nHD=1, nHA=3, tPSA=63.6, Disposed=True, nViolate=0)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 5 - (aMW + alogP + anHD + anHA + atPSA)
    #res
    if detail:
        res = namedtuple('CNS', ['MW', 'logP', 'nHD', 
                                 'nHA', 'tPSA', 'Disposed', 'nViolate'])
        checkres = res(MW, logP, nHD, nHA, tPSA, disposed, violate)
    else:
        res = namedtuple('CNS', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def CheckRespiratory(mol, detail=False):
    """
    Check mol under Respiratory
    
    Ref.:
    -----------
    Ritchie, Timothy J., Christopher N. Luscombe, and Simon JF Macdonald. 
    J Chem Inf Model, 49.4 (2009): 1025-1032.
        
    Rule details:
    -----------
    240<=MW<= 520
    -2.0<=logP<=4.7
    6<=nHB<=12
    51<=tPSA<=135
    3<=nRot<=8
    1<=nRing<=5
    
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    detail: bool(optional, default: False)
        When set to True, function will return more 
        information(the physicochemical properties mentioned in rule)
        else, only return Disposed and nViolated
        
    Return:
    -----------
    checkres: namedtuple
        if detail has been set to True, would be more specific, like:
            Respiratory(MW=180.16, logP=1.31, nHB=4, tPSA=63.6, nRot=3, nRing=1, Disposed=False, nViolate=2)
        else, only return 'Disposed' and 'nViolate', 
        present obey(True) or not the rule and how any requirements has been violated respectively.
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
        disposed = True
    else:
        disposed = False
    #Count the number of violated rules
    violate = 6 - (aMW + alogP + anHB + atPSA + anRot + anRing)
    #res
    if detail:
        res = namedtuple('Respiratory', ['MW', 'logP', 'nHB', 
                                         'tPSA', 'nRot', 'nRing', 
                                         'Disposed', 'nViolate'])
        checkres = res(MW, logP, nHB, tPSA, nRot, nRing, disposed, violate)
    else:
        res = namedtuple('CNS', ['Disposed', 'nViolate'])
        checkres = res(disposed, violate)
    return checkres


def Check_CustomizeRule(mol,prop_kws,closed_interval=True,detail=False):
    """
    You could customize the rule with mostly properties ypu want.
         
    Parameters:
    -----------
    mol: rdkit.Chem.rdchem.Mol
    prop_kws: dict
        the keys of dict are properties you want to check;
        the values should be a tuple or list with two elements,
        present the left- and right-bounded respectively.
    closed_interval: bool, optional(default=True)
        True for using closed interval and False for opened interval
    ---
    Return:
        a namedtuple            
        namedtuple('CheckRes',['Disposed','MatchedAtoms','MatchedNames','Endpoint']), for detail=True;
        and namedtuple('CheckRes',['Disposed','Endpoint']) for False
    """
    allowed = ['MW','nBond','nHet','nRot','nRig',
               'nRing','nHev','logP','MR','nHD',
               'nHA','nHB','AP','logSw','Fsp3','tPSA',
               'AP','HetRatio','MaxRing','nStero', 'HetRatio']    
    properties = molproperty.GetProperties(mol)._asdict()
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
        checkres = namedtuple('CustomizeRule',keys)
    else:
        checkres = namedtuple('CustomizeRule',['nViolate','VioProp'])
        res = [nViolate,VioProp]
    return checkres(*res)
       

#def CheckRule(mol):
#    """
#    """
#    res = [CheckEganRule(mol)[-2],
#           CheckVeberRule(mol)[-2],
#           CheckLipinskiRule(mol)[-2],
#           CheckBeyondRo5(mol)[-2],
#           CheckPfizerRule(mol)[-2],
#           CheckGSKRule(mol)[-2],
#           CheckOralMacrocycles(mol)[-2],
#           CheckOpreaRule(mol)[-2],
#           CheckGhoseRule(mol)[-2],
#           CheckXuRule(mol)[-2],
#           CheckRo4(mol)[-2],
#           CheckRo3(mol)[-2],
#           CheckRo2(mol)[-2],
#           CheckCNS(mol)[-2],
#           CheckRespiratory(mol)[-2]]
#    res = [[False,True][x=='Reject'] for x in res]
#    return res



if __name__ =='__main__':
    
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]','CC(=O)OC1=CC=CC=C1C(=O)O']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
    smiring = ['C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2','C1CCC2C3CC(C4CCCC5CCCCC45)CCC3CCC2C1','C1CCC2(CCC3(CCCCC3)CC2)CC1']
    for index, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        print('Index:{}'.format(index))
        res = CheckRo3(mol, detail=True)
        print(res)
                                  

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   