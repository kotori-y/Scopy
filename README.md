---
@Data: 2019-06-25
@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me
---

# Scopy: a compounds filter for early stages drug discovery.

## what's it

**Scopy** is a Python package, based on RDKit, computing a number of **properties** **of molecule** which present in some **drug likeness rules**, such as Lipinski's rule, Pfizer rule, Beyond Ro5 rule, so that a molecule could be checked under some rules. Besides, Scopy could scan molecule through a predefined toxic fragment in **SMARTS** format to filter the molecule may has some **unexpected endpoint**.

## Main Function

### compute properties

```python
---
up to now(2019.06.24), we have achived followed properties:
    Molcular Weight,MW
    Number of bonds, nBond
    Number of atoms, nAtom
    Number of heteroatoms, nHet
    Number of rotable bonds, nRot
    Number of rigid bonds, nRig
    Number of SSSR, nRing
    Number of heavy atom, nHev
    logP, #we haven't discriminated the AlogP and ClogP
    Molecular refraction, MR
    Number of hydrogen bond donors, nHD
    Number of hydrogen bond acceptors, nHA
    Number of hydrogen bond donors& acceptors, nHB
    Aromatic proportion, AP
    logSw #by the ESOL method
    sp3 hybridized carbons/total carbon count, Fsp3
    TPSA, tpsa
    Number of atoms involved in the biggest system ring, MaxRing
    Number of Sterocenterss, nStero
    HetCarbonRatio, HetRatio
---
Followed should be achieved in the future:
    logD
    difference between clogP and alogP
    Formal total charge of the compound
    Number of Charged Groups
```
this function is realized in module scopy.Druglikeness

```python
>>>from scopy.Druglikeness import CalculateProperty
>>>from rdkit import Chem
```

```python
>>>mol = Chem.MolFromSmiles('C1=CC=CC1')
>>>res = CalculateProperty.GetProperties(mol)
>>>print(res)
Properties(MW=66.1, nBond=5, nAtom=11, nHD=0, nHA=0, nHB=0, nHet=0, nStero=0, nHev=5, nRot=0, nRig=5, nRing=1, logP=1.5, logSw=-1.21, MR=22.9, tPSA=0.0, AP=0.0, HetRatio=0.0, Fsp3=0.2, MaxRing=5)
```

You could also compute each property respectively, like:

```python
>>>mol = Chem.MolFromSmiles('C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2')
>>>res_1 = CalculateProperty.CalculateLogP(mol)
>>>res_2 = CalculateProperty.CalculateNumHAcceptors(mol)
>>>res_3 = CalculateProperty.CalculateNumRing(mol)
>>>print(res_1,res_2,res_3)
5.15 0 4
```

or you could calculate all properties firstly

```python
>>>mol = Chem.MolFromSmiles('C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2')
>>>res = CalculateProperty.GetProperties(mol)
>>>print(res)
>>>print(res.logP,res.nHA,res.nRing)
Properties(MW=228.29, nBond=21, nAtom=30, nHD=0, nHA=0, nHB=0, nHet=0, nStero=0, nHev=18, nRot=0, nRig=21, nRing=4, logP=5.15, logSw=-5.28, MR=78.96, tPSA=0.0, AP=1.0, HetRatio=0.0, Fsp3=0.0, MaxRing=18)
5.15 0 4
```

### Check Druglikeness Rules

```python
---
up to now(2019.06.24), we have achived followed rules:
    Egan Rule
    Veber Rule
    Lipinski Rule
    Beyond Ro5
    Pfizer Rule(3/75 PfizerRule)
    GSK Rule
    Macrocycles
    Oprea Rule
    Ghose Rule
   	Xu Rule
    Ro4
    Ro3
    Ro2
---
Followed should be achieved in the future:
    OpreaTwo Rule
    Kelder Rule
    REOS
    GoldenTriangle rule
    Schneider rule
    DrugLikeOne rule
    DrugLikeTwo rule
    Zinc rule
    CNS rule
    Respiratory rule
```

this function is realized in module scopy.Druglikeness

```python
>>>from scopy.Druglikeness import CheckRule
>>>from rdkit import Chem

>>>mol = Chem.MolFromSmiles('C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2')
>>>res = CheckRule.CheckLipinskiRule(mol)
>>>print(res)
LipinskiRule(Disposed='Accepted', nViolate=1)
```

The filed 'Disposed' is meant molecular state after rule applied. 'Accepted' means obey the rule. attribute 'nViolated' means the number of violated requirement of a rule.

If you want to get specific properties suggested in rule, you could pass the 'True'  to Parameter 'detail'(default: False)

```python
>>>res = DrugLikeness.CheckLipinskiRule(mol,detail=True)
>>>print(res)
LipinskiRule(MW=228.29, logP=5.15, nHD=0, nHA=0, Disposed='Accept', nViolated=1)
```

