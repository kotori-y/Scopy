# -*- coding: utf-8 -*-

#Created on Wed Jul 17 10:15:43 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me


import os
from rdkit import Chem
from rdkit import DataStructs

from .. import ScoConfig

# This module is derived from our previous work
# these are SMARTS patterns corresponding to the PubChem fingerprints
# https://astro.temple.edu/~tua87106/list_fingerprints.pdf
# ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt

class PubChem(object):
    """ This module is derived from our previous work.
    
    These are SMARTS patterns corresponding to the PubChem fingerprints:
        (1) https://astro.temple.edu/~tua87106/list_fingerprints.pdf
        (2) ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
    """
    def __init__(self):
        """Initialization
        
        """
        with open(os.path.join(ScoConfig.PubChemDir, 'pubchem.txt')) as f_obj:
            self.smartsPatts = eval(f_obj.read())
        f_obj.close()    
        self.PubchemKeys = None

    
    def InitKeys(self, keyList, keyDict):
      """ *Internal Use Only*
      
      generates SMARTS patterns for the keys, run once
    
      """
      assert len(keyList) == len(keyDict.keys()), 'length mismatch'
      for key in keyDict.keys():
        patt, count = keyDict[key]
        if patt != '?':
          sma = Chem.MolFromSmarts(patt)
          if not sma:
            print('SMARTS parser error for key #%d: %s' % (key, patt))
          else:
            keyList[key - 1] = sma, count
    
    
    def calcPubChemFingerPart1(self, mol, **kwargs):
      """Calculate PubChem Fingerprints (1-115; 263-881)
      
      :param mol: molecule
      :type mol: rdkit.Chem.rdchem.Mol
      :return: fingerprint
      :rtype: rdkit.DataStructs.cDataStructs.SparseBitVect
      
      """
#      global PubchemKeys
      if self.PubchemKeys is None:
        self.PubchemKeys = [(None, 0)] * len(self.smartsPatts.keys())
    
        self.InitKeys(self.PubchemKeys, self.smartsPatts)
      ctor = kwargs.get('ctor', DataStructs.SparseBitVect)
    
      res = ctor(len(self.PubchemKeys) + 1)
      for i, (patt, count) in enumerate(self.PubchemKeys):
        if patt is not None:
          if count == 0:
            res[i + 1] = mol.HasSubstructMatch(patt)
          else:
            matches = mol.GetSubstructMatches(patt)
            if len(matches) > count:
              res[i + 1] = 1
      return res
    
    
    def func_1(self,mol,bits):
        """ *Internal Use Only*
    
        Calculate PubChem Fingerprints (116-263)
    
        """
        ringSize=[]
        temp={3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        AllRingsAtom = mol.GetRingInfo().AtomRings()
        for ring in AllRingsAtom:
            ringSize.append(len(ring))
            for k,v in temp.items():
                if len(ring) == k:
                    temp[k]+=1
        if temp[3]>=2:
            bits[0]=1;bits[7]=1
        elif temp[3]==1:
            bits[0]=1
        else:
            pass
        if temp[4]>=2:
            bits[14]=1;bits[21]=1
        elif temp[4]==1:
            bits[14]=1
        else:
            pass
        if temp[5]>=5:
            bits[28]=1;bits[35]=1;bits[42]=1;bits[49]=1;bits[56]=1
        elif temp[5]==4:
            bits[28]=1;bits[35]=1;bits[42]=1;bits[49]=1
        elif temp[5]==3:
            bits[28]=1;bits[35]=1;bits[42]=1
        elif temp[5]==2:
            bits[28]=1;bits[35]=1
        elif temp[5]==1:
            bits[28]=1
        else:
            pass
        if temp[6]>=5:
            bits[63]=1;bits[70]=1;bits[77]=1;bits[84]=1;bits[91]=1
        elif temp[6]==4:
            bits[63]=1;bits[70]=1;bits[77]=1;bits[84]=1
        elif temp[6]==3:
            bits[63]=1;bits[70]=1;bits[77]=1
        elif temp[6]==2:
            bits[63]=1;bits[70]=1
        elif temp[6]==1:
            bits[63]=1
        else:
            pass
        if temp[7]>=2:
            bits[98]=1;bits[105]=1
        elif temp[7]==1:
            bits[98]=1
        else:
            pass
        if temp[8]>=2:
            bits[112]=1;bits[119]=1
        elif temp[8]==1:
            bits[112]=1
        else:
            pass
        if temp[9]>=1:
            bits[126]=1;
        else:
            pass
        if temp[10]>=1:
            bits[133]=1;
        else:
            pass
    
        return ringSize,bits
    
    
    def func_2(self,mol,bits):
        """ *Internal Use Only*
    
        saturated or aromatic carbon-only ring
    
        """
        AllRingsBond = mol.GetRingInfo().BondRings()
        ringSize=[]
        temp={3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        for ring in AllRingsBond:
            ######### saturated
            nonsingle = False
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='SINGLE':
                    nonsingle = True
                    break
            if nonsingle == False:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
            ######## aromatic carbon-only     
            aromatic = True
            AllCarb = True
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='AROMATIC':
                    aromatic = False
                    break
            for bondIdx in ring:
                BeginAtom = mol.GetBondWithIdx(bondIdx).GetBeginAtom()
                EndAtom = mol.GetBondWithIdx(bondIdx).GetEndAtom()
                if BeginAtom.GetAtomicNum() != 6 or EndAtom.GetAtomicNum() != 6:
                    AllCarb = False
                    break
            if aromatic == True and AllCarb == True:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
        if temp[3]>=2:
            bits[1]=1;bits[8]=1
        elif temp[3]==1:
            bits[1]=1
        else:
            pass
        if temp[4]>=2:
            bits[15]=1;bits[22]=1
        elif temp[4]==1:
            bits[15]=1
        else:
            pass
        if temp[5]>=5:
            bits[29]=1;bits[36]=1;bits[43]=1;bits[50]=1;bits[57]=1
        elif temp[5]==4:
            bits[29]=1;bits[36]=1;bits[43]=1;bits[50]=1
        elif temp[5]==3:
            bits[29]=1;bits[36]=1;bits[43]=1
        elif temp[5]==2:
            bits[29]=1;bits[36]=1
        elif temp[5]==1:
            bits[29]=1
        else:
            pass
        if temp[6]>=5:
            bits[64]=1;bits[71]=1;bits[78]=1;bits[85]=1;bits[92]=1
        elif temp[6]==4:
            bits[64]=1;bits[71]=1;bits[78]=1;bits[85]=1
        elif temp[6]==3:
            bits[64]=1;bits[71]=1;bits[78]=1
        elif temp[6]==2:
            bits[64]=1;bits[71]=1
        elif temp[6]==1:
            bits[64]=1
        else:
            pass
        if temp[7]>=2:
            bits[99]=1;bits[106]=1
        elif temp[7]==1:
            bits[99]=1
        else:
            pass
        if temp[8]>=2:
            bits[113]=1;bits[120]=1
        elif temp[8]==1:
            bits[113]=1
        else:
            pass
        if temp[9]>=1:
            bits[127]=1;
        else:
            pass
        if temp[10]>=1:
            bits[134]=1;
        else:
            pass
        return ringSize, bits
    
    
    def func_3(self,mol,bits):
        """ *Internal Use Only*
    
        saturated or aromatic nitrogen-containing
    
        """
        AllRingsBond = mol.GetRingInfo().BondRings()
        ringSize=[]
        temp={3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        for ring in AllRingsBond:
            ######### saturated
            nonsingle = False
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='SINGLE':
                    nonsingle = True
                    break
            if nonsingle == False:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
            ######## aromatic nitrogen-containing    
            aromatic = True
            ContainNitro = False
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='AROMATIC':
                    aromatic = False
                    break
            for bondIdx in ring:
                BeginAtom = mol.GetBondWithIdx(bondIdx).GetBeginAtom()
                EndAtom = mol.GetBondWithIdx(bondIdx).GetEndAtom()
                if BeginAtom.GetAtomicNum() == 7 or EndAtom.GetAtomicNum() == 7:
                    ContainNitro = True
                    break
            if aromatic == True and ContainNitro == True:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
        if temp[3]>=2:
            bits[2]=1;bits[9]=1
        elif temp[3]==1:
            bits[2]=1
        else:
            pass
        if temp[4]>=2:
            bits[16]=1;bits[23]=1
        elif temp[4]==1:
            bits[16]=1
        else:
            pass
        if temp[5]>=5:
            bits[30]=1;bits[37]=1;bits[44]=1;bits[51]=1;bits[58]=1
        elif temp[5]==4:
            bits[30]=1;bits[37]=1;bits[44]=1;bits[51]=1
        elif temp[5]==3:
            bits[30]=1;bits[37]=1;bits[44]=1
        elif temp[5]==2:
            bits[30]=1;bits[37]=1
        elif temp[5]==1:
            bits[30]=1
        else:
            pass
        if temp[6]>=5:
            bits[65]=1;bits[72]=1;bits[79]=1;bits[86]=1;bits[93]=1
        elif temp[6]==4:
            bits[65]=1;bits[72]=1;bits[79]=1;bits[86]=1
        elif temp[6]==3:
            bits[65]=1;bits[72]=1;bits[79]=1
        elif temp[6]==2:
            bits[65]=1;bits[72]=1
        elif temp[6]==1:
            bits[65]=1
        else:
            pass
        if temp[7]>=2:
            bits[100]=1;bits[107]=1
        elif temp[7]==1:
            bits[100]=1
        else:
            pass
        if temp[8]>=2:
            bits[114]=1;bits[121]=1
        elif temp[8]==1:
            bits[114]=1
        else:
            pass
        if temp[9]>=1:
            bits[128]=1;
        else:
            pass
        if temp[10]>=1:
            bits[135]=1;
        else:
            pass
        return ringSize, bits
    
    
    def func_4(self,mol,bits):
        """ *Internal Use Only*
    
        saturated or aromatic heteroatom-containing
    
        """
        AllRingsBond = mol.GetRingInfo().BondRings()
        ringSize=[]
        temp={3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        for ring in AllRingsBond:
            ######### saturated
            nonsingle = False
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='SINGLE':
                    nonsingle = True
                    break
            if nonsingle == False:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
            ######## aromatic heteroatom-containing
            aromatic = True
            heteroatom = False
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='AROMATIC':
                    aromatic = False
                    break
            for bondIdx in ring:
                BeginAtom = mol.GetBondWithIdx(bondIdx).GetBeginAtom()
                EndAtom = mol.GetBondWithIdx(bondIdx).GetEndAtom()
                if BeginAtom.GetAtomicNum() not in [1,6] or EndAtom.GetAtomicNum() not in [1,6]:
                    heteroatom = True
                    break
            if aromatic == True and heteroatom == True:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
        if temp[3]>=2:
            bits[3]=1;bits[10]=1
        elif temp[3]==1:
            bits[3]=1
        else:
            pass
        if temp[4]>=2:
            bits[17]=1;bits[24]=1
        elif temp[4]==1:
            bits[17]=1
        else:
            pass
        if temp[5]>=5:
            bits[31]=1;bits[38]=1;bits[45]=1;bits[52]=1;bits[59]=1
        elif temp[5]==4:
            bits[31]=1;bits[38]=1;bits[45]=1;bits[52]=1
        elif temp[5]==3:
            bits[31]=1;bits[38]=1;bits[45]=1
        elif temp[5]==2:
            bits[31]=1;bits[38]=1
        elif temp[5]==1:
            bits[31]=1
        else:
            pass
        if temp[6]>=5:
            bits[66]=1;bits[73]=1;bits[80]=1;bits[87]=1;bits[94]=1
        elif temp[6]==4:
            bits[66]=1;bits[73]=1;bits[80]=1;bits[87]=1
        elif temp[6]==3:
            bits[66]=1;bits[73]=1;bits[80]=1
        elif temp[6]==2:
            bits[66]=1;bits[73]=1
        elif temp[6]==1:
            bits[66]=1
        else:
            pass
        if temp[7]>=2:
            bits[101]=1;bits[108]=1
        elif temp[7]==1:
            bits[101]=1
        else:
            pass
        if temp[8]>=2:
            bits[115]=1;bits[122]=1
        elif temp[8]==1:
            bits[115]=1
        else:
            pass
        if temp[9]>=1:
            bits[129]=1;
        else:
            pass
        if temp[10]>=1:
            bits[136]=1;
        else:
            pass
        return ringSize,bits
    
    
    def func_5(self,mol,bits):
        """ *Internal Use Only*
    
        unsaturated non-aromatic carbon-only
    
        """
        ringSize=[]
        AllRingsBond = mol.GetRingInfo().BondRings()
        temp={3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        for ring in AllRingsBond:
            unsaturated = False
            nonaromatic = True
            Allcarb = True
            ######### unsaturated
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='SINGLE':
                    unsaturated = True
                    break
            ######## non-aromatic
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name=='AROMATIC':
                    nonaromatic = False
                    break
            ######## allcarb
            for bondIdx in ring:
                BeginAtom = mol.GetBondWithIdx(bondIdx).GetBeginAtom()
                EndAtom = mol.GetBondWithIdx(bondIdx).GetEndAtom()
                if BeginAtom.GetAtomicNum() != 6 or EndAtom.GetAtomicNum() != 6:
                    Allcarb = False
                    break
            if unsaturated == True and nonaromatic == True and Allcarb == True:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
        if temp[3]>=2:
            bits[4]=1;bits[11]=1
        elif temp[3]==1:
            bits[4]=1
        else:
            pass
        if temp[4]>=2:
            bits[18]=1;bits[25]=1
        elif temp[4]==1:
            bits[18]=1
        else:
            pass
        if temp[5]>=5:
            bits[32]=1;bits[39]=1;bits[46]=1;bits[53]=1;bits[60]=1
        elif temp[5]==4:
            bits[32]=1;bits[39]=1;bits[46]=1;bits[53]=1
        elif temp[5]==3:
            bits[32]=1;bits[39]=1;bits[46]=1
        elif temp[5]==2:
            bits[32]=1;bits[39]=1
        elif temp[5]==1:
            bits[32]=1
        else:
            pass
        if temp[6]>=5:
            bits[67]=1;bits[74]=1;bits[81]=1;bits[88]=1;bits[95]=1
        elif temp[6]==4:
            bits[67]=1;bits[74]=1;bits[81]=1;bits[88]=1
        elif temp[6]==3:
            bits[67]=1;bits[74]=1;bits[81]=1
        elif temp[6]==2:
            bits[67]=1;bits[74]=1
        elif temp[6]==1:
            bits[67]=1
        else:
            pass
        if temp[7]>=2:
            bits[102]=1;bits[109]=1
        elif temp[7]==1:
            bits[102]=1
        else:
            pass
        if temp[8]>=2:
            bits[116]=1;bits[123]=1
        elif temp[8]==1:
            bits[116]=1
        else:
            pass
        if temp[9]>=1:
            bits[130]=1;
        else:
            pass
        if temp[10]>=1:
            bits[137]=1;
        else:
            pass
        return ringSize,bits
    
    
    def func_6(self,mol,bits):
        """ *Internal Use Only*
    
        unsaturated non-aromatic nitrogen-containing
    
        """
        ringSize=[]
        AllRingsBond = mol.GetRingInfo().BondRings()
        temp={3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        for ring in AllRingsBond:
            unsaturated = False
            nonaromatic = True
            ContainNitro = False
            ######### unsaturated
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='SINGLE':
                    unsaturated = True
                    break
            ######## non-aromatic
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name=='AROMATIC':
                    nonaromatic = False
                    break
            ######## nitrogen-containing
            for bondIdx in ring:
                BeginAtom = mol.GetBondWithIdx(bondIdx).GetBeginAtom()
                EndAtom = mol.GetBondWithIdx(bondIdx).GetEndAtom()
                if BeginAtom.GetAtomicNum() == 7 or EndAtom.GetAtomicNum() == 7:
                    ContainNitro = True
                    break
            if unsaturated == True and nonaromatic == True and ContainNitro== True:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
        if temp[3]>=2:
            bits[5]=1;bits[12]=1
        elif temp[3]==1:
            bits[5]=1
        else:
            pass
        if temp[4]>=2:
            bits[19]=1;bits[26]=1
        elif temp[4]==1:
            bits[19]=1
        else:
            pass
        if temp[5]>=5:
            bits[33]=1;bits[40]=1;bits[47]=1;bits[54]=1;bits[61]=1
        elif temp[5]==4:
            bits[33]=1;bits[40]=1;bits[47]=1;bits[54]=1
        elif temp[5]==3:
            bits[33]=1;bits[40]=1;bits[47]=1
        elif temp[5]==2:
            bits[33]=1;bits[40]=1
        elif temp[5]==1:
            bits[33]=1
        else:
            pass
        if temp[6]>=5:
            bits[68]=1;bits[75]=1;bits[82]=1;bits[89]=1;bits[96]=1
        elif temp[6]==4:
            bits[68]=1;bits[75]=1;bits[82]=1;bits[89]=1
        elif temp[6]==3:
            bits[68]=1;bits[75]=1;bits[82]=1
        elif temp[6]==2:
            bits[68]=1;bits[75]=1
        elif temp[6]==1:
            bits[68]=1
        else:
            pass
        if temp[7]>=2:
            bits[103]=1;bits[110]=1
        elif temp[7]==1:
            bits[103]=1
        else:
            pass
        if temp[8]>=2:
            bits[117]=1;bits[124]=1
        elif temp[8]==1:
            bits[117]=1
        else:
            pass
        if temp[9]>=1:
            bits[131]=1;
        else:
            pass
        if temp[10]>=1:
            bits[138]=1;
        else:
            pass
        return ringSize,bits
    
    
    def func_7(self,mol,bits):
        """ *Internal Use Only*
    
        unsaturated non-aromatic heteroatom-containing
    
        """
        ringSize=[]
        AllRingsBond = mol.GetRingInfo().BondRings()
        temp={3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        for ring in AllRingsBond:
            unsaturated = False
            nonaromatic = True
            heteroatom = False
            ######### unsaturated
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='SINGLE':
                    unsaturated = True
                    break
            ######## non-aromatic
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name=='AROMATIC':
                    nonaromatic = False
                    break
            ######## heteroatom-containing
            for bondIdx in ring:
                BeginAtom = mol.GetBondWithIdx(bondIdx).GetBeginAtom()
                EndAtom = mol.GetBondWithIdx(bondIdx).GetEndAtom()
                if BeginAtom.GetAtomicNum() not in [1,6] or EndAtom.GetAtomicNum() not in [1,6]:
                    heteroatom = True
                    break
            if unsaturated == True and nonaromatic == True and heteroatom == True:
                ringSize.append(len(ring))
                for k,v in temp.items():
                    if len(ring) == k:
                        temp[k]+=1
        if temp[3]>=2:
            bits[6]=1;bits[13]=1
        elif temp[3]==1:
            bits[6]=1
        else:
            pass
        if temp[4]>=2:
            bits[20]=1;bits[27]=1
        elif temp[4]==1:
            bits[20]=1
        else:
            pass
        if temp[5]>=5:
            bits[34]=1;bits[41]=1;bits[48]=1;bits[55]=1;bits[62]=1
        elif temp[5]==4:
            bits[34]=1;bits[41]=1;bits[48]=1;bits[55]=1
        elif temp[5]==3:
            bits[34]=1;bits[41]=1;bits[48]=1
        elif temp[5]==2:
            bits[34]=1;bits[41]=1
        elif temp[5]==1:
            bits[34]=1
        else:
            pass
        if temp[6]>=5:
            bits[69]=1;bits[76]=1;bits[83]=1;bits[90]=1;bits[97]=1
        elif temp[6]==4:
            bits[69]=1;bits[76]=1;bits[83]=1;bits[90]=1
        elif temp[6]==3:
            bits[69]=1;bits[76]=1;bits[83]=1
        elif temp[6]==2:
            bits[69]=1;bits[76]=1
        elif temp[6]==1:
            bits[69]=1
        else:
            pass
        if temp[7]>=2:
            bits[104]=1;bits[111]=1
        elif temp[7]==1:
            bits[104]=1
        else:
            pass
        if temp[8]>=2:
            bits[118]=1;bits[125]=1
        elif temp[8]==1:
            bits[118]=1
        else:
            pass
        if temp[9]>=1:
            bits[132]=1;
        else:
            pass
        if temp[10]>=1:
            bits[139]=1;
        else:
            pass
        return ringSize,bits
    
    
    def func_8(self,mol,bits):
        """ *Internal Use Only*
    
        aromatic rings or hetero-aromatic rings
    
        """
        AllRingsBond = mol.GetRingInfo().BondRings()
        temp={'aromatic':0,'heteroatom':0}
        for ring in AllRingsBond:
            aromatic = True
            heteroatom = False
            for bondIdx in ring:
                if mol.GetBondWithIdx(bondIdx).GetBondType().name!='AROMATIC':
                    aromatic = False
                    break
            if aromatic==True:
                temp['aromatic']+=1
            for bondIdx in ring:
                BeginAtom = mol.GetBondWithIdx(bondIdx).GetBeginAtom()
                EndAtom = mol.GetBondWithIdx(bondIdx).GetEndAtom()
                if BeginAtom.GetAtomicNum() not in [1,6] or EndAtom.GetAtomicNum() not in [1,6]:
                    heteroatom = True
                    break
            if heteroatom==True:
                temp['heteroatom']+=1
        if temp['aromatic']>=4:
            bits[140]=1;bits[142]=1;bits[144]=1;bits[146]=1
        elif temp['aromatic']==3:
            bits[140]=1;bits[142]=1;bits[144]=1
        elif temp['aromatic']==2:
            bits[140]=1;bits[142]=1
        elif temp['aromatic']==1:
            bits[140]=1
        else:
            pass
        if temp['aromatic']>=4 and temp['heteroatom']>=4:
            bits[141]=1;bits[143]=1;bits[145]=1;bits[147]=1
        elif temp['aromatic']==3 and temp['heteroatom']==3:
            bits[141]=1;bits[143]=1;bits[145]=1
        elif temp['aromatic']==2 and temp['heteroatom']==2:
            bits[141]=1;bits[143]=1
        elif temp['aromatic']==1 and temp['heteroatom']==1:
            bits[141]=1
        else:
            pass
        return bits
    
    
    def calcPubChemFingerPart2(self,mol):# 116-263
        """ *Internal Use Only*
    
        Calculate PubChem Fingerprints （116-263)
    
        """
        bits=[0]*148
        bits=self.func_1(mol,bits)[1]
        bits=self.func_2(mol,bits)[1]
        bits=self.func_3(mol,bits)[1]
        bits=self.func_4(mol,bits)[1]
        bits=self.func_5(mol,bits)[1]
        bits=self.func_6(mol,bits)[1]
        bits=self.func_7(mol,bits)[1]
        bits=self.func_8(mol,bits)
    
        return bits
    
    
    def CalculatePubChem(self,mol):
        """Calculate PubChem Fingerprints
        
        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol
        :return: fingerprint
        :rtype: list
        
        """
        AllBits=[0]*881
        res1=list(self.calcPubChemFingerPart1(mol).ToBitString())
        for index, item in enumerate(res1[1:116]):
            if item == '1':
                AllBits[index] = 1
        for index2, item2 in enumerate(res1[116:734]):
            if item2 == '1':
                AllBits[index2+115+148] = 1
        res2=self.calcPubChemFingerPart2(mol)
        for index3, item3 in enumerate(res2):
            if item3==1:
                AllBits[index3+115]=1
        return AllBits


if __name__ == '__main__':
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

    mol = Chem.MolFromSmiles(smis[3])
    fp = PubChem()
    res = fp.CalculatePubChem(mol)
#    fps = np.array(fps)
    print(res)