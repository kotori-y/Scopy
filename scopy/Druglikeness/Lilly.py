# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 21:59:42 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

"""
##################################################################################################
This script is aim at achieving check molecule with Lilly’s rulers under Python3 environment;

-Ref:
    Bruns, Robert F., and Ian A. Watson.
    Journal of medicinal chemistry 55.22 (2012): 9763-9772.
    
version == 1.0
##################################################################################################
"""

from rdkit.Chem import AllChem as Chem
from rdkit.Chem.rdchem import BondType
import time
"""
########################################################################################
The following atomic attributes are recognized within a query file
########################################################################################
"""

def _Check_atomic_number(atom,AtomicNum):
    """
    **Internal used only**

    -element specification    
    """    
    num = atom.GetAtomicNum()
    return num in AtomicNum


def _Check_atomic_symbol(atom,AtomSymbol):
    """
    **Internal used only** 
    
    -element specification
    """  
    pass


def _Check_ncon(atom,con):
    """
    **Internal used only** 
    
    number of explicit connections
    """ 
    num = len(atom.GetBonds())
    return num == con


def _Check_min_ncon(atom,con):
    """
    **Internal used only** 
    
    number of explicit connections
    """ 
    num = len(atom.GetBonds())
    return num >= con


def _Check_ncon2(atom):
    """
    **Internal used only**

    number explicit connections within 2 bonds    
    """ 
    num = len(atom.GetBonds())
    return num <=2


def _Check_nbonds(atom,bondnum):
    """
    **Internal used only**

    number of bonds (1 == single, 2 == double…)    
    """ 
    num = 0
    for bond in atom.GetBonds():
        btype = bond.GetBondType()
        if btype == BondType.SINGLE:
            num += 1
        elif btype == BondType.TRIPLE:
            num += 3
        elif btype == BondType.AROMATIC or btype == BondType.DOUBLE:
            num += 2
    return num == bondnum         
    

def _Check_min_nbonds(atom,minnum):
    """
    **Internal used only**

    number of bonds (1 == single, 2 == double…)    
    """ 
    num = 0
    for bond in atom.GetBonds():
        btype = bond.GetBondType()
        if btype == BondType.SINGLE:
            num += 1
        elif btype == BondType.TRIPLE:
            num += 3
        elif btype == BondType.AROMATIC or btype == BondType.DOUBLE:
            num += 2
    return num >= minnum       


def _Check_formal_charge(atom,chargenum):
    """
    **Internal used only** 
    
    formal charge specification
    """ 
    num = atom.GetFormalCharge()
    return num == chargenum
    
    
def _Check_nrings(mol,atom,numsssr):
    """
    **Internal used only**  
    
    number of SSSR rings
    """ 
    num = 0
    idx = atom.GetIdx()
    info = mol.GetRingInfo()
    ringatoms = info.AtomRings()
    for a in ringatoms:
        if a == idx:
            num +=1
    return num == numsssr


def _Check_ring_bond_count(atom,numcount):
    """
    **Internal used only**

    number of connections via ring bonds    
    """ 
    pass
    

def _Check_ring_size(atom,size):
    """
    **Internal used only** 
    
    size of a ring containing atom
    """ 
    return atom.IsInRingSize(size)
 
    
def _Check_aromatic_ring_size(atom,size):
    """
    **Internal used only** 
    
    size of an aromatic ring containing atom
    """ 
    if atom.GetIsAromatic():
        return atom.IsInRingSize(size)
    else:
        if size != 0:
            return False
        else:
            return True


def _Check_aliphatic_ring_size(atom,size):
    """
    **Internal used only**  
    
    size of an aliphatic ring containing atom
    """ 
    if not atom.GetIsAromatic():
        return atom.IsInRingSize(size)
    else:
        return False
    

def _Check_hcount(atom,count):
    """
    **Internal used only**

    number of explicit + implicit hydrogens    
    """ 
    explih = atom.GetNumExplicitHs()
    implih = atom.GetNumImplicitHs()
    return explih+implih


def _Check_attached_heteroatom_count(atom,numheter):
    """
    **Internal used only**

    number of heteroatoms attached  
    """ 
    neighbors = atom.GetNeighbors()
    neighbors = [atom.GetAtomicNum() for atom in neighbors\
                 if atom.GetAtomicNum() not in [1,6]]
    return len(neighbors) == numheter
    
    
def _Check_lone_pair(atom):
    """
    **Internal used only**

    number of lone pairs  
    """
    pass
    

def _Check_unsaturation(atom,numcount):
    """
    **Internal used only**

    nbonds – ncon  
    """
    nbonds = 0
    for bond in atom.GetBonds():
        btype = bond.GetBondType()
        if btype == BondType.SINGLE:
            nbonds += 1
        elif btype == BondType.TRIPLE:
            nbonds += 3
        elif btype == BondType.AROMATIC or btype == BondType.DOUBLE:
            nbonds += 2
    ncon = len(atom.GetBonds())
    return (nbonds - ncon) == numcount


def _Check_daylight_x(atom,numx):
    """
    **Internal used only**

    ncon + implicit hydrogens
    """
    ncon = len(atom.GetBonds())
    implih = atom.GetNumImplicitHs()
    return (ncon + implih) == numx


def _Check_isotope(atom,numiso):
    """
    **Internal used only**

    isotopic label – can be any positive number
    """
    pass


def _Check_aryl(atom,arnum):
    """
    **Internal used only**

    adjacent to an aromatic ring
    """
    neighbors = atom.GetNeighbors()
    neighbors = [neighbor for neighbor in neighbors\
                 if neighbor.GetIsAromatic()]
    return len(neighbors) == arnum


def _Check_vinyl(atom,vinum):
    """
    **Internal used only**

    adjacent to a non-aromatic, unsaturated atom.
    """
    num = 0
    neighbors = atom.GetNeighbors()
    neighbors = [neighbor for neighbor in neighbors\
                 if not neighbor.GetIsAromatic()]
    for neighbor in neighbors:
        if neighbor.IsInRing():
            num += 1
        else:
            for bond in neighbor.GetBonds():
                btype = bond.GetBondType()
                if btype == BondType.DOUBLE or btype == BondType.TRIPLE:
                    num += 1
                    break
    return num == vinum
    
    
def _Check_fused_system_size(mol,atom,fusedsize):
    """
    **Internal used only**

    size of fused system containing atom
    """
    num = 0
    sta = []
    flag = 1
    info = mol.GetRingInfo()
    ringatoms = list(info.AtomRings())
    tar = [ringatoms.pop(idx) for idx,x in enumerate(ringatoms) if atom.GetIdx() in x]
    while flag == 1:
        for idx,ring in enumerate(ringatoms):
            if len(set(ring)&set(tar)) >= 2:
                num += 1
                sta.append(ringatoms.pop(idx))
                break
            else:
                flag = 0
    
    """
    I give up temporarily
    """
    pass

            
def _Check_heteroatoms_in_ring(mol,atom,heternum):
    """
    **Internal used only**

    in ring containing this many heteroatoms
    """  
    num = 0
    info = mol.GetRingInfo()
    ringatoms = list(info.AtomRings())
    tar = [x.GetAtomicNum() for x in ringatoms if atom.GetIdx() in x]
    for item in tar[0]:
        if item != 6 and item != 1:
            num += 1
    return num == heternum


def _Check_bond_exist(atom,btype):
    """
    **Internal used only**

    Checking specific bond whether exist or not
    1 == SINGLE
    2 == DOUBLE
    3 == TRIPLE
    4 == AROMATIC    
    """  
    bonds = atom.GetBonds()
    if btype == 1:
        res = [bond for bond in bonds if bond.GetBondType() == BondType.SINGLE]
    elif btype == 2:
        res = [bond for bond in bonds if bond.GetBondType() == BondType.DOUBLE]
    elif btype == 3:
        res = [bond for bond in bonds if bond.GetBondType() == BondType.TRIPLE]
    elif btype == 4:
        res = [bond for bond in bonds if bond.GetBondType() == BondType.AROMATIC]
    return len(res)>=1
        
"""
########################################################################################
The following Ring System Specification are recognized within a query file
########################################################################################
"""

def _Get_nrings(mol):
    """
    **Internal used only**

    number of SSSR rings in the ring system
    """  
    return Chem.GetSSSR(mol) 
    
    
def _Get_aromatic_ring_count(ringinfo):
    """
    **Internal used only**

    number of aromatic rings
    """
    num = 0
    for item in ringinfo.AtomRings():
        if mol.GetAtoms()[item[0]].GetIsAromatic():
            num += 1
    return num

def _Get_ring_size(mol):
    """
    **Internal used only**

    number atoms in the ring
    """
    info = mol.GetRingInfo()
    sizes = [len(x) for x in (info.AtomRings())]
    return sizes



##########################################Followed are queries##########################################


def Check_3_valent_halogen(atoms):    
#    Followed is query file in original Ruby script:
#        
#        (0 Query
#          (A I Version 2)
#          (A C Comment "3_valent_halogen")
#          (0 Query_Atom
#            (A I atomic_number (9 17 35 53))
#            (A I nbonds 3)
#          )
#        )        
    for atom in atoms:
        if _Check_atomic_number(atom,AtomicNum=[9,17,35,53]):
            if _Check_nbonds(atom,bondnum=3):
                return 'Reject'
    return 'Accept'
            

def Check_3_valent_iodine(atoms):
#     Followed is query file in original Ruby script:
#         
#        (0 Query
#          (A I Version 2)
#          (A C Comment "3_valent_iodine")
#          (0 Query_Atom
#            (A I atomic_number 53)
#            (A I nbonds 3)
#          )
#        )
    for atom in atoms:
        if _Check_atomic_number(atom,AtomicNum=[53]):
            if _Check_nbonds(atom,bondnum=3):
                return 'Reject'
    return 'Accept'
          

def Check_4_valent_sulphur_2_connections(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A I Version 2)
#          (A C Comment "4_valent_sulphur_2_connections")
#          (A C smarts "[SD2v4]")
#        )
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[SD2v4]')):
        return 'Reject'
    else:
        return 'Accept'

    
def Check_6_membered_aromatic_sulfur(mol):
#    Followed is query file in original Ruby script:
#        (0 Querynnric_value 50)
#          (A C Comment "6_membered_aromatic_sulfur")
#          (A C smarts "s1aaaaa1")
#        )
    """
    Demerits
    """
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts('s1aaaaa1'))) >= 2:
        return 'Reject'
    else:
        return 'Accept'
    

def Check_8_aminoquinoline(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A C Comment "8_aminoquinoline")
#          (A I min_nrings 2)
#          (A I Version 2)
#          (0 ring_system_specifier
#            (A I nrings 2)
#            (A I aromatic_ring_count 2)
#          )
#          (A C smarts "[nD2H0]1c2c(-[ND2H])cccc2ccc1")
#        )
    """
    should be revised
    """
    numring = _Get_nrings(mol)
    if numring < 2:
        return 'Accept'
    else:
        if _Get_aromatic_ring_count(mol.GetRingInfo()) < 2:
            return 'Accept'
        else:
            if mol.HasSubstructMatch(Chem.MolFromSmarts('[nD2H0]1c2c(-[ND2H])cccc2ccc1')):
                return 'Reject'
            else:
                return 'Accept'
        

def Check_8_hydroxyquinoline(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A I Version 2)
#          (A C Comment "8_hydroxyquinoline")
#          (A C smarts "[OD1]c1cccc2cccnc12")
#        )
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[OD1]c1cccc2cccnc12')):
        return 'Reject'
    else:
        return 'Accept' 


def Check_9_aminoacridine(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A C Comment "8_aminoquinoline")
#          (A I min_nrings 2)
#          (A I Version 2)
#          (0 ring_system_specifier
#            (A I nrings 2)
#            (A I aromatic_ring_count 2)
#          )
#          (A C smarts "[nD2H0]1c2c(-[ND2H])cccc2ccc1")
#        )
    """
    should be revised
    """
    numring = _Get_nrings(mol)
    if numring < 3:
        return 'Accept'
    else:
        if _Get_aromatic_ring_count(mol.GetRingInfo()) < 3:
            return 'Accept'
        else:
            if mol.HasSubstructMatch(Chem.MolFromSmarts('N-c1c2c([nH0]c3c1cccc3)cccc2')):
                return 'Reject'
            else:
                return 'Accept'


def Check_acetal_1_in_ring(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A I Version 2)
#          (A D numeric_value 30)
#          (A C Comment "acetal_1_in_ring")
#          (A C smarts "C[O,S;R1]C[O,S;R0]C")
#        )
    """
    Demerits
    
    This may should be revised in future
    """
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts('C[O,S;R1]C[O,S;R0]C'))) >= 4:
        return 'Reject'
    else:
        return 'Accept'
            

def Check_acetal_acyclic(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A C Comment "acetal_acyclic")
#          (A C Comment "acetal_acyclic")
#          (A I Version 2)
#          (0 Query_Atom
#            (A I atomic_number (7 8 16))
#            (A I ncon 2)
#            (A I nrings 0)
#          )
#          (1 Query_Atom
#            (A I atomic_number 6)
#            (A I unsaturation 0)
#            (A I nrings 0)
#            (A I single_bond 0)
#          )
#          (2 Query_Atom
#            (A I atomic_number (7 8 16))
#            (A I ncon 2)
#            (A I nrings 0)
#            (A I single_bond 1)
#          )
#        )
    atoms = mol.GetAtoms()
    for atom in atoms:
        #0 Query_Atom
        if _Check_atomic_number(atom,AtomicNum=[7,8,16]) and\
        _Check_ncon(atom,con=2) and\
        _Check_nrings(mol,atom,numsssr=0):
            #1 Query_Atom
            for atom in atoms:
                if _Check_atomic_number(atom,AtomicNum=[6]) and\
                _Check_unsaturation(atom,0) and\
                _Check_nrings(mol,atom,0) and\
                len([bond for bond in atom.GetBonds() if bond.GetBondType() == BondType.SINGLE])>=1:
                    bondidx_0 = [bond.GetIdx() for bond in atom.GetBonds() if bond.GetBondType() == BondType.SINGLE]
                    #2 Query_Atom
                    for atom in atoms:                            
                        if _Check_atomic_number(atom,AtomicNum=[7,8,16]) and\
                        _Check_ncon(atom,2) and\
                        _Check_nrings(mol,atom,0) and\
                        len([bond for bond in atom.GetBonds() if bond.GetBondType() == BondType.SINGLE])>=1:                     
                            bondidx_1 = [bond.GetIdx() for bond in atom.GetBonds() if bond.GetBondType() == BondType.SINGLE]
                            if len(set(bondidx_1)^set(bondidx_0)) > 0:
                                return 'Reject'
                    return 'Accept'
            return'Accept'
    return 'Accept'
            
                            
def Check_acetal_both_in_ring(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A I Version 2)
#          (A D numeric_value 30)
#          (A C Comment "acetal_both_in_ring")
#          (A C smarts "C[O,S;R1]C[O,S;R1]C")
#        )
    """
    Demerits
    """
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts('C[O,S;R1]C[O,S;R1]C'))) >= 4:
        return 'Reject'
    else:
        return 'Accept'
             
        
def Check_acetate_ester(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A C Comment "acetate_ester")
#          (A I Version 2)
#          (A I min_hits_needed 3)
#          (0 Query_Atom
#            (A I atomic_number 8)
#            (A I ncon 2)
#            (A I attached_heteroatom_count 0)
#          )
#          (1 Query_Atom
#            (A I atomic_number 6)
#            (A I ncon 3)
#            (A I nbonds 4)
#            (A I attached_heteroatom_count 2)
#            (A I single_bond 0)
#          )
#          (2 Query_Atom
#            (A I atomic_number 8)
#            (A I double_bond 1)
#          )
#          (3 Query_Atom
#            (A I atomic_number 6)
#            (A I ncon 1)
#            (A I unsaturation 0)
#            (A I single_bond 1)
#          )
#        )
    count = 0
    atoms = mol.GetAtoms()
    for atom in atoms:
        #0 Query
        if _Check_atomic_number(atom,[8]) and\
        _Check_ncon(atom,2) and\
        _Check_attached_heteroatom_count(atom,0):
            count += 1
            if count >=3:
                count = 0
                for atom in atoms:
                    #1 Query
                    if _Check_atomic_number(atom,[6]) and\
                    _Check_ncon(atom,3) and\
                    _Check_nbonds(atom,4) and\
                    _Check_attached_heteroatom_count(atom,2):
                        count += 1
                        if count >=3:
                            count = 0
                            for atom in atoms:
                                #2 Query
                                if _Check_atomic_number(atom,[8]) and\
                                len([bond for bond in atom.GetBonds() if bond.GetBondType() == BondType.DOUBLE])>=1:
                                    count += 1
                                    if count >= 3:
                                        count = 0
                                        for atom in atoms:
                                            if _Check_atomic_number(atom,[6]) and\
                                            _Check_ncon(atom,1) and\
                                            _Check_unsaturation(atom,0):
                                                count += 1
                                                if count >= 3:
                                                    return 'Reject'
                                        return 'Accept'
                            return 'Accept'
                return 'Accept'
    return 'Accept'


def Check_acetylene(mol):
#    Followed is query file in original Ruby script:
#        (3 Query
#          (A C Comment "acetylene")
#          (A I Version 2)
#          (A I unique_embeddings_only 1)
#          (A D numeric_value 50)
#          (0 Query_Atom
#            (A I atomic_number 6)
#            (A I min_nbonds 3)
#            (A I nrings 0)
#          )
#          (1 Query_Atom
#            (A I atomic_number 6)
#            (A I min_ncon 2)
#            (A I nrings 0)
#            (A I triple_bond 0)
#          )
#        )
    """
    Demerits
    """
    count = 0
    atoms = mol.GetAtoms()
    for atom in atoms:
        #0 Query
        if _Check_atomic_number(atom,[6]) and\
        _Check_min_nbonds(atom,3) and\
        _Check_nrings(mol,atom,0):
            count += 1
            if count >= 2:
                count = 0
                for atom in atoms:
                    #1 Query
                    if _Check_atomic_number(atom,[6]) and\
                    _Check_min_ncon(atom,2) and\
                    _Check_nrings(mol,atom,0) and\
                    len([bond for bond in atom.GetBonds() if bond.GetBondType() == BondType.TRIPLE])>=1:
                        count += 1
                        if count >=2:
                            return 'Reject'
                return 'Accept'
    return 'Accept'
        
    
def Check_acetylene_heteroatom(mol):
#    Followed is query file in original Ruby script:
#        (4 Query
#          (A C Comment "acetylene_heteroatom")
#          (A I Version 2)
#          (0 Query_Atom
#            (A I atomic_number 6)
#            (A I ncon 2)
#            (A I nbonds 4)
#            (A I attached_heteroatom_count 1)
#          )
#          (1 Query_Atom
#            (A I atomic_number 6)
#            (A I triple_bond 0)
#          )
#        )
    atoms = mol.GetAtoms()
    for atom in atoms:
        #0 Query
        if _Check_atomic_number(atom,[6]) and\
        _Check_ncon(atom,2) and\
        _Check_nbonds(atom,4) and\
        _Check_attached_heteroatom_count(atom,1):
            for atom in atoms:
                if _Check_atomic_number(atom,[6]) and\
                _Check_bond_exist(atom,3):
                    return 'Reject'
            return 'Accept'
    return 'Accept'
            

def Check_acid_halide(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A I Version 2)
#          (A C Comment "acid_halide")
#          (A I one_embedding_per_start_atom 1)
#          (0 Query_Atom
#            (A I atomic_number (9 17 35 53))
#            (A I ncon 1)
#          )
#          (1 Query_Atom
#            (A I atomic_number 6)
#            (A I ncon 3)
#            (A I nbonds 4)
#            (A I single_bond 0)
#          )
#          (2 Query_Atom
#            (A I atomic_number 8)
#            (A I ncon 1)
#            (A I double_bond 1)
#          )
#        )
    atoms = mol.GetAtoms()
    for atom in atoms:
        #0 Query
        if _Check_atomic_number(atom,[9,17,35,53]) and\
        _Check_ncon(atom,1):
            for atom in atoms:
                #1 Query
                if _Check_atomic_number(atom,[6]) and\
                _Check_ncon(atom,3) and\
                _Check_nbonds(atom,4) and\
                _Check_bond_exist(atom,1):
                    for atom in atoms:
                        #2 Query
                        if _Check_atomic_number(atom,[8]) and\
                        _Check_ncon(atom,1) and\
                        _Check_bond_exist(atom,2):
                            return 'Reject'
                    return 'Accept'
            return 'Accept'
    return 'Accept'


def Check_activated_ester(mol):
#    Followed is query file in original Ruby script:
#        (0 Query
#          (A I Version 2)
#          (A C Comment "activated_ester")
#          (A I unique_embeddings_only 1)
#          (0 Query_Atom
#            (A I atomic_number 8)
#            (A I ncon 1)
#            (A I nbonds 2)
#          )
#          (1 Query_Atom
#            (A I atomic_number (6 16 15))
#            (A I min_ncon 3)
#            (A I nrings 0)
#            (A I double_bond 0)
#          )
#          (2 Query_Atom
#            (A I atomic_number 8)
#            (A I ncon 2)
#            (A I single_bond 1)
#          )
#          (3 Query_Atom
#            (A I atomic_number (7 16 15))
#            (A I single_bond 2)
#          )
#        )
    atoms = mol.GetAtoms()
    for atom in atoms:
        #0 Query
        if _Check_atomic_number(atom,[8]) and\
        _Check_ncon(atom,1) and\
        _Check_nbonds(atom,2):
            for atom in atoms:
                #1 Query
                if _Check_atomic_number(atom,[6,16,15]) and\
                _Check_min_ncon(atom,3) and\
                _Check_nrings(mol,atom,0) and\
                _Check_bond_exist(atom,2):
                    for atom in atoms:
                        #2 Query
                        if _Check_atomic_number(atom,[8]) and\
                        _Check_ncon(atom,2) and\
                        _Check_bond_exist(atom,1):
                            bondidx_0 = [bond.GetIdx() for bond in atom.GetBonds() if bond.GetBondType() == BondType.SINGLE]
                            for atom in atoms:
                                #3 Query
                                if _Check_atomic_number(atom,[7,16,15]) and\
                                _Check_bond_exist(atom,1):
                                    bondidx_1 = [bond.GetIdx() for bond in atom.GetBonds() if bond.GetBondType() == BondType.SINGLE]
                                    if len(set(bondidx_1)^set(bondidx_0)) >= 1:
                                        return 'Reject'                                
                            return 'Accept'
                    return 'Accept'
                return 'Accept'
    return  'Accept'
    
    
def Check_activated_phthalimide(mol):
#    Followed is query file in original Ruby script:
#        (12 Query
#          (A C Comment "activated_phthalimide")
#          (A I Version 2)
#          (A I min_nrings 2)
#          (A I embeddings_do_not_overlap 1)
#          (14 Ring_System_Specifier
#            (A I nrings 2)
#            (A I ring_size (5 6))
#          )
#          (0 Query_Atom
#            (A I atomic_number 8)
#            (A I ncon 1)
#            (A I nbonds 2)
#            (A I nrings 0)
#          )
#          (1 Query_Atom
#            (A I atomic_number 6)
#            (A I ncon 3)
#            (A I nbonds 4)
#            (A I nrings 1)
#            (A I ring_size 5)
#            (A I aromatic 0)
#            (A I ring_id 1)
#            (A I double_bond 0)
#          )
#          (2 Query_Atom
#            (A I atomic_number 7)
#            (A I nrings 1)
#            (A I ring_size 5)
#            (A I aromatic 0)
#            (A I ring_id 1)
#            (A I single_bond 1)
#          )
#          (3 Query_Atom
#            (A I atomic_number 6)
#            (A I ncon 3)
#            (A I nbonds 4)
#            (A I nrings 1)
#            (A I ring_size 5)
#            (A I aromatic 0)
#            (A I ring_id 1)
#            (A I single_bond 2)
#          )
#          (4 Query_Atom
#            (A I atomic_number 8)
#            (A I ncon 1)
#            (A I nbonds 2)
#            (A I nrings 0)
#            (A I aromatic 0)
#            (A I double_bond 3)
#          )
#          (5 Query_Atom
#            (A I atomic_number 6)
#            (A I ncon 3)
#            (A I nbonds 4)
#            (A I nrings 2)
#            (A I ring_size (5 6))
#            (A I aromatic 1)
#            (A I single_bond 3)
#          )
#          (6 Query_Atom
#            (A I atomic_number 6)
#            (A I min_ncon 2)
#            (A I nrings 1)
#            (A I ring_size 6)
#            (A I aromatic 1)
#            (A I aromatic_bond 5)
#          )
#          (7 Query_Atom
#            (A I atomic_number 6)
#            (A I min_ncon 2)
#            (A I nrings 1)
#            (A I ring_size 6)
#            (A I aromatic 1)
#            (A I aromatic_bond 6)
#          )
#          (8 Query_Atom
#            (A I atomic_number 6)
#            (A I min_ncon 2)
#            (A I nrings 1)
#            (A I ring_size 6)
#            (A I aromatic 1)
#            (A I aromatic_bond 7)
#          )
#          (9 Query_Atom
#            (A I atomic_number 6)
#            (A I min_ncon 2)
#            (A I nrings 1)
#            (A I ring_size 6)
#            (A I aromatic 1)
#            (A I aromatic_bond 8)
#          )
#          (10 Query_Atom
#            (A I atomic_number 6)
#            (A I ncon 3)
#            (A I nbonds 4)
#            (A I nrings 2)
#            (A I ring_size (5 6))
#            (A I aromatic 1)
#            (A I aromatic_bond 9)
#            (A I aromatic_bond 5)
#            (A I single_bond 1)
#          )
#          (19 Environment
#            (A I single_bond (6 7 8 9))
#            (A C smarts "[F,Cl,Br,I]")
#            (A C smiles "N(=O)=O")
#            (A C smarts "C=O")
#          )
#        )
    if _Get_nrings(mol)>=2 and\
    5 in _Get_ring_size(mol) and\
    6 in _Get_ring_size(mol):
        atoms = mol.GetAtoms()
        for atom in atoms:
            #0 Query_Atom
            if _Check_atomic_number(atom,[8]) and\
            _Check_ncon(atom,1) and\
            _Check_nbonds(atom,2) and\
            _Check_nrings(atom,0):
                for atom in atoms:
                    #1 Query_Atom
                    if _Check_atomic_number(atom,6) and\
                    _Check_ncon(atom,3) and\
                    _Check_nbonds(atom,4) and\
                    _Check_nrings(atom,1) and\
                    _Check_ring_size(atom,5) and\
                    _Check_bond_exist(atom,2) and\
                    atom.GetIsAromatic():
                        info = mol.GetRingInfo()
                        ringatoms_a = [set(x) for x in info.AtomRings() if atom.GetIdx() in x]
                        bondidx_a_d = [bond.GetIdx() for bond in atom.GetBonds() if bond.GetBondType() == BondType.DOUBLE]
                        for atom in atoms:
                            #2 Query_Atom
                            if _Check_atomic_number(atom,[7]) and\
                            _Check_nrings(atom,1) and\
                            _Check_ring_size(atom,5) and\
                            _Check_bond_exist(atom,1) and\
                            atom.GetIsAromatic():
                               ringatoms_b = [set(x) for x in info.AtomRings() if atom.GetIdx() in x]
                               for item in ringatoms_b:
                                   if item in ringatoms_a:
                                      bondidx_a_s = [bond.GetIdx() for bond in atom.GetBonds() if bond.GetBondType() == BondType.SINGLE]
                                      for atom in atoms:
                                          #3 Query_Atom
                                          pass
                        return 'Accept'
                return 'Accept'
                    
        return 'Accept'            
    else:
        return 'Accept'
      
"""
To be continued...
"""    
    
    
    
    
    
    

if '__main__' == __name__:
    smis = [
            'C1(=CC=CC=C1)I=O', #3_valent_halogen
            'C1(=CC=CC=C1)I=O', #3_valent_iodine
            'C1(=O)C(=S=NO1)C', #4_valent_sulphur_2_connections
            'C1(=NSN=C(Cl)C1=O)N1CCN(C2=NSN=C(Cl)C2=O)CC1', #6_membered_aromatic_sulfur
            'C12=NC=CC=C1C=CC=C2NN', #8_aminoquinoline
            'C12=CC=CN=C1C(=CC=C2)O', #8_hydroxyquinoline
            'C12=CC=CC=C1C(=C1C=CC=CC1=N2)N', #9_aminoacridine
            'C1(=NCCS1)SC', #acetal_1_in_ring
            'C(Br)[C@@H](OC)OC', #acetal_acyclic
            'C1(=NC#N)SCCS1', #acetal_both_in_ring
            'C([C@@H](COC(C)=O)OC(C)=O)OC(C)=O', #acetate_ester
            'N(CC#C)CC#C', #acetylene
            'C(#CCC)OCC', #acetylene_heteroatom
            'C(C)(C)(C)C(=O)Cl', #acid_halide
            'C(=NOC(C)=O)(C)C', #activated_ester
            '', #activated_phthalimide
            'C1=CC=CC(/C=C/C#C)=C1'
            ]
    
    smi = smis[-1]
    mol = Chem.MolFromSmiles(smi)
    #atoms = mol.GetAtoms()
#    for atom in mol.GetAtoms():
#        idx = atom.GetIdx()
#        res = _Check_fused_system_size(mol,atom,1)
#        print(idx,res)
    start = time.process_time()      
    print(Check_activated_phthalimide(mol))
    end = time.process_time()    
    print(end-start)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        