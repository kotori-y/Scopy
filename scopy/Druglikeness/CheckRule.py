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
    This module implments properties obtained from module CalculateProperty
    ---
    We have collected followed rules:

    Egan Rule     0<=tPSA<=132; -1<=logP<=6
    Veber Rule   nRot<= 10; tPSA<= 140; nHB<= 12
    LipinskiRule  MW<=500; logP<=5, nHD<=5, nHA <=10
    BeyondRo5   Mw<=1000; -2<=logP<=10; nHD<=6, nHA<=15; tPSA<=250; nRot<=20
    Pfizer Rule     logP>3; PSA<75
    GSK Rule        MW<=400; logP<=4
    OralMacrocycles     MW<1000; logP<10; nHD<5; PSA<250
    Oprea Rule      0<=nHD<=2; 2<= nHA<=9; 2<=nRot<=8; 1<=nRing<=4
    Ghose Rule      -0.4<logP<5.6; 160<MW<480; 40<MR<130; 20<nAtom<70
    Xu Rule     nHD<=5; nHA <= 10; 3 <= rot <= 35; 1 <= nring <= 7; 10 <= nhev <= 50
    Ro4 Rule    MW<=400; logP<=4; nHD<=4; NHA<=8; PSA<=120
    Ro3 Rule    MW<=300; -3<=logP<=3; nHD<=3; nHA<=6; PSA<=60
    Ro2 Rule    MW<=200; logP<=2; nHD<=2; nHA<=4

    ---
    Followed should be achieved in the future:

    OpreaTwo Rule       mw <= 450; -3.5<=clogP<=4.5; -4<=logD<=4; nring<= 4; no. of nonterminal single bonds<=10; nHD<=5; nHA<=8
    Kelder Rule
    REOS Rule
    GoldenTriangle
    Schneider Rule
    DrugLikeOne
    DrugLikeTwo
    Zinc
    Cns
    Respiratory
    """

from scopy.Druglikeness import CalculateProperty
from collections import namedtuple
from rdkit import Chem

def CheckEganRule(mol,detail):
    """
    #################################################################
    Bad or Good oral biovailability rule
    
    -Ref.:
        Egan, William J., Kenneth M. Merz, and John J. Baldwin. 
        J Med Chem, 43.21 (2000): 3867-3877.
        
    -Rule details:
        0 <= tPSA <= 132
        -1 <= logP <=6
    #################################################################
    """
    #Calculate required properties of Lipinskin's rule
    tPSA = CalculateProperty.CalculateTPSA(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    
    #Determine whether the molecular match each rule
    atPSA = (0<=tPSA<=132)
    alogP = (-1<=logP<=6)   
    #Give the advice
    if atPSA&alogP:
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violsted rules
    violate = 2 - (atPSA+alogP)    
    #res
    if detail:
        res = namedtuple('EganRuler',['tPSA','logP','Disposed','Violate'])
        checkres = res(tPSA,logP,disposed,violate)
    else:
        res = namedtuple('EganRuler',['Disposed','Violate'])
        checkres = res(disposed,violate)
    return checkres


def CheckVeberRule(mol,detail=False):
    """
    #################################################################
    Bad or Good oral biovailability rule
    
    -Ref.:
        Veber, Daniel F., et al.
        Journal of medicinal chemistry 45.12 (2002): 2615-2623.
        
    -Rule details:
        nRot <= 10
        TPSA <= 140
        the number of H-Bonds Donors and H-Bonds Acceptors <= 12
    #################################################################
    """
    #Calculate required properties of Lipinskin's rule
    nRot = CalculateProperty.CalculateNumRotatableBonds(mol)
    tPSA = CalculateProperty.CalculateTPSA(mol)
    nHB = CalculateProperty.CalculateNumHyBond(mol)   
    #Determine whether the molecular match each rule
    anRot = (nRot<=10)
    atPSA = (tPSA<=140)
    anHB = (nHB<=12)    
    #Give the advice
    if (anRot&atPSA&anHB):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 3 - (anRot+atPSA+anHB)
    #res
    if detail:
        res = namedtuple('VeberRuler',['nRot','tPSA','nHB','Disposed','Violate'])
        checkres = res(nRot,tPSA,nHB,disposed,violate)
    else:
        res = namedtuple('VeberRuler',['Disposed','Violate'])
        checkres = res(disposed,violate)    
    return checkres


def CheckLipinskiRule(mol,detail=False):
    """
    #################################################################
    Check molecular under Lipinski's rule
    
    -Ref.:
        Lipinski, Christopher A., et al.
        Advanced drug delivery reviews 23.1-3 (1997): 3-25.
    
    -Rule details:
        molecular weight <= 500(Da)
        LogP <= 5
        number of hydrogen bond donors <= 5
        number of hydrogen Bond accptors <= 10
    #################################################################
    """
    #Calculate required properties of Lipinskin's rule
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    nHA = CalculateProperty.CalculateNumHAcceptors(mol)   
    #Determine whether the molecular match each rule
    aw = (MW <= 500)
    alogp = (logP <= 5)
    anhd = (nHD <= 5)
    anha = (nHA <= 10)
    #Count the number of matched rules
    nviolate = 4 - (aw+alogp+anhd+anha)
    #Give the disposed
    if nviolate >= 2:
        disposed = 'Reject'
    else:
        disposed = 'Accept'
        
    if detail:
        res = namedtuple('LipinskiRule',['MW','logP','nHD','nHA','Disposed','nViolate'])
        checkres = res(MW,logP,nHD,nHA,disposed,nviolate)
    else:
        res = namedtuple('LipinskiRule',['Disposed','nViolate'])
        checkres = res(disposed,nviolate)    
    return checkres   
    

def CheckBeyondRo5(mol,detail=False):
    """
    #################################################################
    Check molecular under beyond Ro5
    
    -Ref.:
        Doak, Bradley C., et al.
        journal of medicinal chemistry 59.6 (2015): 2312-2327.
        
    -Rule details:
        MW <= 1000
        -2 <= logP <= 10
        nHD <= 6
        nHA <= 15
        PSA <= 250
        nRot <= 20        
    #################################################################
    """
    #Calculate required properties of Lipinskin's rule
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    nHA = CalculateProperty.CalculateNumHAcceptors(mol)
    PSA = CalculateProperty.CalculateTPSA(mol)
    nRot = CalculateProperty.CalculateNumRotatableBonds(mol)
    
    #Determine whether the molecular match each rule
    aMW = (MW<=1000)
    alogP = (-2<=logP<=10)
    anHD = (nHD<=6)
    anHA = (nHA<=15)
    aPSA = (PSA<=250)
    aRot = (nRot<=20)
    
    #Check whether mol pass the whole rules
    if (aMW&alogP&anHD&anHA&aRot&aPSA):
        disposed = 'Accept'
    else:
        disposed = 'Reject'   
    #Count the number of violated rules
    violate = 6 - (aMW+alogP+anHD+anHA+aRot+aPSA)    
    if detail:
        res = namedtuple('bRo5Ruler',['MW','logP','nHD','nHA','PSA','nRot','Disposed','Violate'])
        checkres = res(MW,logP,nHD,nHA,PSA,nRot,disposed,violate)
    else:
        res = namedtuple('bRo5Ruler',['Disposed','Violate'])
        checkres = res(disposed,violate)        
    return checkres    


def CheckPfizerRule(mol,detail=False):
    """
    #################################################################
    Check molecular under Rfizer Rule(3/75 Rule)
    
    -Ref.:
        Hughes, Jason D., et al. 
        Bioorganic & medicinal chemistry letters 18.17 (2008): 4872-4875.
        
    -Rule details:
        ---Logp > 3
        ---PSA < 75
    #################################################################
    """
    #Calculate required properties
    logP = CalculateProperty.CalculateLogP(mol)
    PSA = CalculateProperty.CalculateTPSA(mol)
    
    #Determine whether the molecular match each rule
    alogP = (logP > 3)
    aPSA = (PSA < 75)    
    #Give the advice
    if alogP&aPSA:
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 2 - (alogP+aPSA)
    #res
    if detail:
        res = namedtuple('PfizerRule',['logP','PSA','Disposed','Violate'])
        checkres = res(logP,PSA,disposed,violate)
    else:
        res = namedtuple('PfizerRule',['Disposed','Violate'])
        checkres = res(disposed,violate)    
    return checkres
    

def CheckGSKRule(mol,detail=False):
    """
    #################################################################
    Check molecular under GSK rule(4/400 Rule)
    
    -Ref.:
        Gleeson, M. Paul.
        Journal of medicinal chemistry 51.4 (2008): 817-834.
        
    -Rule details:
        molecular weight <= 400
        LogP <= 4
    #################################################################
    """
    #Calculate required properties
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)   
    #Determine whether the molecular match each rule    
    aMW = (MW <= 400)
    alogP = (logP <= 4)    
    #Give the advice
    if aMW&alogP:
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 2 - (aMW+alogP)
    #res
    if detail:
        res = namedtuple('GSKRuler',['MW','logP','Disposed','Violate'])
        checkres = res(MW,logP,disposed,violate)
    else:
        res = namedtuple('GSKRuler',['Disposed','Violate'])
        checkres = res(disposed,violate)
    return checkres


def OralMacrocycles(mol,detail=False):
    """
    #################################################################
    Check molecular under oral macrocycles rules
    -Ref.:
        Giordanetto, Fabrizio, and Jan Kihlberg.
        Journal of medicinal chemistry 57.2 (2013): 278-295.        

    -Rule details:
        molecular weight < 1000(Da)
        LogP < 10
        number of hydrogen bond donors < 5
        PSA < 250
    #################################################################
    """
    #Calculate required properties of Lipinskin's rule
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    tPSA = CalculateProperty.CalculateTPSA(mol)   
    #Determine whether the molecular match each rule
    aMW = (MW < 1000)
    alogP = (logP < 10)
    anHD = (nHD < 5)
    aPSA = (tPSA < 250)  
    #Give the advice
    if (aMW&alogP&anHD&aPSA):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violate rules
    violate = 4 - (aMW+alogP+anHD+aPSA)
    #res
    if detail:
        res = namedtuple('OralMacrocycles',['MW','logP','nHD','tPSA','Disposed','Violate'])
        checkres = res(MW,logP,nHD,tPSA,disposed,violate)
    else:
        res = namedtuple('OralMacrocycles',['MW','logP','nHD','tPSA','Disposed','Violate'])
        checkres = res(MW,logP,nHD,tPSA,disposed,violate) 
    return checkres


def CheckOpreaRule(mol,detail=False):
    """
    #################################################################
    -Ref.:
        Oprea, Tudor I.
        Journal of computer-aided molecular design 14.3 (2000): 251-264.
        
    -Rule details:
        0 <= nHD <= 2
        2 <= nHA <= 9
        2 <= nRot <= 8
        1 <= nRing <= 4
    #################################################################
    """
    #Calculate required properties of Oprea
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    nHA = CalculateProperty.CalculateNumHAcceptors(mol)
    nRot = CalculateProperty.CalculateNumRotatableBonds(mol)
    nRing = CalculateProperty.CalculateNumRing(mol)    
    #Determine whether the molecular match each rule
    anHD = (0<=nHD<=2)
    anHA = (2<=nHA<=9)
    anRot = (2<=nRot<=8)
    anRing = (1<=nRing<=4)   
    #Give the advice
    if (anHD&anHA&anRot&anRing):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 4 - (anHD+anHA+anRot+anRing)
    #res   
    if detail:
        res = namedtuple('OpreaRules',['nHD','nHA','nRot','nRing','Disposed','Violate'])
        checkres = res(nHD,nHA,nRot,nRing,disposed,violate)
    else:
        res = namedtuple('OpreaRules',['Disposed','Violate'])
        checkres = res(disposed,violate)
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


def CheckGhoseRule(mol,detail=False):
    """
    #################################################################
    Check molecular under Ghose rule
    
    -Ref.:
        Ghose, Arup K., Vellarkad N. Viswanadhan, and John J. Wendoloski. 
        Journal of combinatorial chemistry 1.1 (1999): 55-68.
    
    -Rules details:
        -0.4<logP<5.6
        160<MW<480
        40<MR<130
        20<nAtom<70
    #################################################################
    """
    #Calculate required properties of Oprea
    logP = CalculateProperty.CalculateLogP(mol)
    MW = CalculateProperty.CalculateMolWeight(mol)
    MR = CalculateProperty.CalculateMolMR(mol)
    nAtom = CalculateProperty.CalculateNumAtoms(mol)    
    #Determine whether the molecular match each rule
    alogP = (-0.4<logP<5.6)
    aMW = (160<MW<480)
    aMR = (40<MR<130)
    anAtom = (20<nAtom<70)
    #Give the advice
    if (alogP&aMW&aMR&anAtom):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 4 - (alogP+aMW+aMR+anAtom)
    #Res
    if detail:
        res = namedtuple('GhoseRuler',['logP','MW','MR','nAtom','Disposed','Violate'])
        checkres = res(logP,MW,MR,nAtom,disposed,violate)
    else:
        res = namedtuple('GhoseRuler',['Disposed','Violate'])
        checkres = res(disposed,violate)
    return checkres


def CheckKelderRule(mol):
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


def CheckREOS(mol):
    """
    #################################################################
    Check molecular under REOS program
    
    -Ref.:
        Walters, W. Patrick, and Mark Namchuk.
        Nat Rev Drug Discov, 2.4 (2003): 259.
        
    -Rule details:
        200 <= mw <= 500
        -5 <= logP <= 5
        nhd <= 5
        nha <= 10
        nrot <= 8
        tpsa <= 150
        -4 <= totalchar <= 4
    #################################################################
    """
    #Calculate required properties of REOS
    mw = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nhd = CalculateProperty.CalculateNumHDonors(mol)
    nha = CalculateProperty.CalculateNumHAcceptors(mol)
    nrot = CalculateProperty.CalculateNumRotatableBonds(mol)
    tpsa = CalculateProperty.CalculateTPSA(mol)
    totalchar = CalculateProperty.CalculateTotalCharge(mol)
    
    #Determine whether the molecular match each rule
    aw = ((mw>=200)&(mw<=500))
    alogP = ((logP>=-5)&(logP<=5))
    anhd = (nhd<=5)
    anha = (nha<=10)
    anrot = ((nrot>=2)&(nrot<=8))
    atpsa = (tpsa<=150)
    atotalchar = ((totalchar>=-4)&(totalchar<=4))

#    #Check whether mol pass the whole rules
#    temp1 = (aw&alogP&anhd&anha&anrot&atpsa&atotalchar)
#    #Count the number of matched rules
#    temp2 = (aw+alogP+anhd+anha+anrot+atpsa+atotalchar)
#    
#    res = [temp1,temp2]
#    
#    return res
    
    
    pass

    return None


def CheckGoldenTriangle(mol):
    """
    #################################################################
    Check molecular under 'Golden Triangle'
    
    NOTICE: We should use LogD!!
    
    -Ref.:
        Johnson, Ted W., Klaus R. Dress, and Martin Edwards.
        Bioorg Med Chem Lett, 19.19 (2009): 5560-5564.
        
    -Rule details:
        200 <= MW <= 500
        -2 <= LogD <= 5
    #################################################################
    """
    #Calculate required properties of REOS
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    
    #Determine whether the molecular match each rule
    aw = (200<=MW<=500)
    alogp = (-2<=logD<=5)
    
    #Check whether mol pass the whole rules
    temp1 = (aw&alogp)
    #Count the number of matched rules
    temp2 = (aw+alogp)
    
    res = [temp1,temp2]
    
    pass


def CheckXuRule(mol,detail):
    """
    #################################################################
    Check molecular under Xu's rule
    
    -Ref.:
        
    -Rule details:
        nhd <= 5
        nha <= 10
        3 <= rot <= 35
        1 <= nring <= 7
        10 <= nhev <= 50
    
    #################################################################
    """
    #Calculate required properties of Xu's rule
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    nHA = CalculateProperty.CalculateNumHAcceptors(mol)
    nRot = CalculateProperty.CalculateNumRotatableBonds(mol)
    nRing = CalculateProperty.CalculateNumRing(mol)
    nHev = CalculateProperty.CalculateHeavyAtomNumber(mol)   
    #Determine whether the molecular match each rule
    anHD = (nHD<=5)
    anHA = (nHA<=10)
    anRot = (3<=nRot<=35)
    anRing = (1<=nRing<=7)
    anHev = (10<=nHev<=50)
    #Give the advice
    if (anHD&anHA&anRot&anRing&anHev):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 5 - (anHD+anHA+anRot+anRing+anHev)
    #res
    if detail:
        res = namedtuple('XuRule',['nHD','nHA','nRot','nRing','nHev','Disposed','Violate'])
        checkres = res(nHD,nHA,nRot,nRing,nHev,disposed,violate)
    else:
        res = namedtuple('XuRule',['Disposed','Violate'])
        checkres = res(disposed,violate)
    return checkres


def CheckSchneiderRule(mol):
    """
    #################################################################
    Check  molecular under Schneider rule
    
    -Ref.:
        Schneider, Nadine, et al.
        J Chem Inf Model, 48.3 (2008): 613-628.
        
    -Rule details:
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


def CheckRo4(mol,detail=False):
    """
    #################################################################
    
    -Ref.:
        
    -Rule:
        MW <= 400
        Logp <= 4
        NHD <= 4
        NHA <= 8
        PSA <= 120
    #################################################################
    """
    #Calculate required properties of REOS
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    nHA = CalculateProperty.CalculateNumHAcceptors(mol)
    tPSA = CalculateProperty.CalculateTPSA(mol)    
    #Determine whether the molecular match each rule
    aMW = (MW<=400)
    alogP = (logP<=4)
    anHD = (nHD<=4)
    anHA = (nHA<=8)
    atPSA = (tPSA<=120)   
    #Give the advice
    if (aMW&alogP&anHD&anHA&atPSA):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 5 - (anHD+anHA+aMW+alogP+atPSA)
    #res
    if detail:
        res = namedtuple('Ro4',['MW','logP','nHD','nHA','tPSA','Disposed','Violate'])
        checkres = res(MW,logP,nHD,nHA,tPSA,disposed,violate)
    else:
        res = namedtuple('Ro4',['Disposed','Violate'])
        checkres = res(disposed,violate)
    return checkres


def CheckRo3(mol,detail=False):
    """
    #################################################################
    Check molecular under Ro3
    
    -Ref.:
        Congreve, Miles, et al.
        Drug discovery today 19.8 (2003): 876-877.
        
    -Rule details:
        MW <= 300
        -3 <= logP <= 3
        NHD <= 3
        NHA <= 6
        PSA <= 60
    #################################################################
    """
    #Calculate required properties of REOS
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    nHA = CalculateProperty.CalculateNumHAcceptors(mol)
    tPSA = CalculateProperty.CalculateTPSA(mol)
    nRot = CalculateProperty.CalculateNumRotatableBonds(mol)    
    #Determine whether the molecular match each rule
    aMW = (MW<=300)
    alogP = (-3<=logP<=3)
    anHD = (nHD<=3)
    anHA = (nHA<=6)
    atPSA = (tPSA<=60)
    anRot = (nRot<=3)
    #Give the advice
    if (aMW&alogP&anHD&anHA&atPSA&anRot):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 6 - (aMW+alogP+anHD+anHA+atPSA+anRot)
    #res
    if detail:
        res = namedtuple('Ro3',['MW','logP','nHD','nHA','tPSA','nRot','Disposed','Violate'])
        checkres = res(MW,logP,nHD,nHA,tPSA,nRot,disposed,violate)
    else:
        res = namedtuple('Ro3',['Disposed','Violate'])
        checkres = res(disposed,violate)    
    return checkres


def CheckRo2(mol,detail):
    """
    #################################################################
    Check molecular under RO2
    
    -Ref.:
        Goldberg, Frederick W., et al.
        Drug Discovery Today 20.1 (2015): 11-17.
        
    -Rule details:
        MW <= 200
        Logp <= 2
        NHD <= 2
        NHA <= 4
    #################################################################    
    """
    #Calculate required properties of REOS
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    nHA = CalculateProperty.CalculateNumHAcceptors(mol)    
    #Determine whether the molecular match each rule
    aMW = (MW<=200)
    alogP = (logP<=2)
    anHD = (nHD<=2)
    anHA = (nHA<=4)    
    #Give the advice
    if (aMW&alogP&anHD&anHA):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 4 - (aMW+alogP+anHD+anHA)
    #res
    if detail:
        res = namedtuple('Ro2',['MW','logP','nHD','nHA','Disposed','Violate'])
        checkres = res(MW,logP,nHD,nHA,disposed,violate)
    else:
        res = namedtuple('Ro2',['Disposed','Violate'])
        checkres = res(disposed,violate)
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


def CheckCNS(mol,detail=False):
    """
    #################################################################
    Check mol under CNS
    
    -Ref.:
        Jeffrey, Phil, and Scott Summerfield.
        Neurobiol Dis, 37.1 (2010): 33-37.
        
    -Rule details:
        135 <= weight <= 582
        -0.2 <= logP <= 6.1
        nha <= 5
        nhd <= 3
        3 <= TPSA <= 118
    #################################################################
    """
    #Calculate required properties of REOS
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nHD = CalculateProperty.CalculateNumHDonors(mol)
    nHA = CalculateProperty.CalculateNumHAcceptors(mol)
    tPSA = CalculateProperty.CalculateTPSA(mol)
    #Determine whether the molecular match each rule
    aMW = (135<=MW<=582)
    alogP = (-0.2<=logP<=6.1)
    anHD = (nHD<=3)
    anHA = (nHA<=5)  
    atPSA = (3<=tPSA<=118)
    #Give the advice
    if (aMW&alogP&anHD&anHA&atPSA):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 5 - (aMW+alogP+anHD+anHA+atPSA)
    #res
    if detail:
        res = namedtuple('CNS',['MW','logP','nHD','nHA','tPSA','Disposed','Violate'])
        checkres = res(MW,logP,nHD,nHA,tPSA,disposed,violate)
    else:
        res = namedtuple('CNS',['Disposed','Violate'])
        checkres = res(disposed,violate)
    return checkres


def CheckRespiratory(mol,detail=False):
    """
    #################################################################
    Check mol under Respiratory
    
    -Ref.:
        Ritchie, Timothy J., Christopher N. Luscombe, and Simon JF Macdonald. 
        J Chem Inf Model, 49.4 (2009): 1025-1032.
        
    -Rule details:
        240<=MW<= 520
        -2.0<=logP<=4.7
        6<=nHB<=12
        51<=tPSA<=135
        3<=nRot<=8
        1<=nRing<=5
    #################################################################
    """
    #Calculate required properties of REOS
    MW = CalculateProperty.CalculateMolWeight(mol)
    logP = CalculateProperty.CalculateLogP(mol)
    nHB = CalculateProperty.CalculateNumHyBond(mol)
    tPSA = CalculateProperty.CalculateTPSA(mol)
    nRot = CalculateProperty.CalculateNumRotatableBonds(mol)
    nRing = CalculateProperty.CalculateNumRing(mol)
    #Determine whether the molecular match each rule
    aMW = (240<=MW<=520)
    alogP = (-0.2<=logP<=4.7)
    anHB = (6<=nHB<=12)  
    atPSA = (51<=tPSA<=135)
    anRot = (3<=nRot<=8)
    anRing = (1<=nRing<=5)
    #Give the advice
    if (aMW&alogP&anHB&atPSA&anRot&anRing):
        disposed = 'Accept'
    else:
        disposed = 'Reject'
    #Count the number of violated rules
    violate = 6 - (aMW+alogP+anHB+atPSA+anRot+anRing)
    #res
    if detail:
        res = namedtuple('Respiratory',['MW','logP','nHB','tPSA','nRot','nRing','Disposed','Violate'])
        checkres = res(MW,logP,nHB,tPSA,nRot,nRing,disposed,violate)
    else:
        res = namedtuple('CNS',['Disposed','Violate'])
        checkres = res(disposed,violate)
    return checkres


def Check_CustomizeRule(mol,prop_kws,closed_interval=True,detail=False):
    """
    You could customize the rule you want.
    
    ---     
    Parameters:
    >>> mol: rdkit.Chem.rdchem.Mol
    >>> prop_kws: dict, the keys of dict are properties you want to check;
                        the values should be a tuple or list with two elements,
                        present the left- and right-bounded respectively.
    >>> closed_interval: bool, optional(default=True), 
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
    properties = CalculateProperty.GetProperties(mol)._asdict()
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
       

if __name__ =='__main__':
    
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]','CC(=O)OC1=CC=CC=C1C(=O)O']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
    smiring = ['C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2','C1CCC2C3CC(C4CCCC5CCCCC45)CCC3CCC2C1','C1CCC2(CCC3(CCCCC3)CC2)CC1']
    for index, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        print('Index:{}'.format(index))
        res = CheckRespiratory(mol,detail=True)
#        print(res)
#        res = Check_CustomizeRule(mol,
#                                  prop_kws={'MW':(100,500),
#                                            'nRot':(2,3),
#                                            'nHD':(None,5),
#                                            'nHB':(1,None)},
#                                  detail=True,
#                                  closed_interval=False)
        print(res)
                                  

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   