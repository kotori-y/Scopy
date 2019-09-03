---
@Data: 2019-06-25
@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China，
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me
---

# Scopy: a compounds filter for early stages drug discovery.

## what's its

**Scopy** is a Python package, based on RDKit, computing various **properties of molecule** which present in some **drug likeness rules**, such as Lipinski's rule, Pfizer rule, Beyond Ro5 rule, so that a molecule could be checked under some rules. Besides, Scopy could scan molecule through a predefined toxic fragments in **SMARTS** format to filter the molecule which may has some **unexpected endpoint**.

## Main Function

### Compute properties

```python
"""
This moudle is used to calculated properise that contained in our collectded rules
    ---
    up to now(2019.07.22), we have achived followed properties:
        Molcular Weight >>> MW
        Number of bonds >>> nBond
        Number of atoms >>> nAtom
        Number of heteroatoms >>> nHet
        Number of heavy atom >>> nHev
        Number of rotable bonds >>> nRot
        Number of rigid bonds >>> nRig
        Number of SSSR >>> nRing
        logP >>> logP
        logD >>> logD
        logSw >>> logSw
        Acid or Base >>> Acid/Base
        pKa >>> pKa
        QED >>> qed
        Molecular refraction >>> MR
        Number of hydrogen bond donors >>> nHD
        Number of hydrogen bond acceptors >>> nHA
        Number of hydrogen bond donors& acceptors >>> nHB
        Aromatic proportion >>> AP
        sp3 hybridized carbons/total carbon count >>> Fsp3
        TPSA >>> TPSA
        Number of atoms involved in the biggest system ring >>> MaxRing
        Number of Sterocenterss >>> nStero
        HetCarbonRatio >>> HetRatio
        synthetic accessibility score >>> SAscore
        natural product-likeness score >>> NPscore
        Number of single bonds >>> nSingle
        Number of double bobds >>> nDoudle
        Number of triple bonds >>> nTriple
        Volume of mol >>> Vol
        Density >>> Dense
        MolFCharge >>> fChar
        Number of Carbon atoms
        Number of Boron atoms
        Number of Chlorin atoms
        Number of Bromine atoms
        Number of Iodine atoms
        Number of Phosphor atoms
        Number of Sulfur atoms
        Number of Oxygen atoms
        Number of Nitrogen atoms       
    ---
    Followed should be achieved in the future:
        logD
        difference between clogP and alogP
        Formal total charge of the compound
        Number of Charged Groups
"""
```
this function is realized in module scopy.Druglikeness

```python
>>> from scopy.Druglikeness import CalculateProperty
>>> from rdkit import Chem
```

```python
>>> mol = Chem.MolFromSmiles('C1=CC=CC1')
>>> res = CalculateProperty.GetProperties(mol)
>>> print(res)
Properties(MW=66.1, nBond=5, nAtom=11, nCarbon=5, nHD=0, nHA=0, nHB=0, nHet=0, nStero=0, nHev=5, nRot=0, nRig=5, nRing=1, logP=1.5, logSw=-1.21, MR=22.9, tPSA=0.0, AP=0.0, HetRatio=0.0, Fsp3=0.2, MaxRing=5, QED=0.4, SAscore=3.31, NPscore=2.16)
```

You could also compute each property respectively, like:

```python
>>> mol = Chem.MolFromSmiles('C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2')
>>> res_1 = CalculateProperty.CalculateLogP(mol)
>>> res_2 = CalculateProperty.CalculateNumHAcceptors(mol)
>>> res_3 = CalculateProperty.CalculateNumRing(mol)
>>> print(res_1,res_2,res_3)
5.15 0 4
```

or you could calculate all properties firstly

```python
>>> mol = Chem.MolFromSmiles('C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2')
>>> res = CalculateProperty.GetProperties(mol)
>>> print(res)
5.15 0 4
>>> print(res.logP,res.nHA,res.nRing)
Properties(MW=228.29, nBond=21, nAtom=30, nCarbon=18, nHD=0, nHA=0, nHB=0, nHet=0, nStero=0, nHev=18, nRot=0, nRig=21, nRing=4, logP=5.15, logSw=-5.28, MR=78.96, tPSA=0.0, AP=1.0, HetRatio=0.0, Fsp3=0.0, MaxRing=18, QED=0.37, SAscore=1.35, NPscore=-0.14)
```

### Check Drug-likeness Rules

```python
"""
---
up to now(2019.07.09), we have achived followed rules:
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
    CNS rule
    Respiratory rule`
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
"""
```

this function is realized in module scopy.Druglikeness

```python
>>> from scopy.Druglikeness import CheckRule
>>> from rdkit import Chem

>>> mol = Chem.MolFromSmiles('C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2')
>>> res = CheckRule.CheckLipinskiRule(mol)
>>> print(res)
LipinskiRule(Disposed='Accepted', nViolate=1)
```

The filed 'Disposed' is meant molecular state after rule applied. 'Accepted' means obey the rule. attribute 'nViolated' means the number of violated requirement of a rule.

If you want to get specific properties suggested in rule, you could pass the <code>True</code>  to Parameter 'detail'(default: False)

```python
>>> res = CheckRule.CheckLipinskiRule(mol,detail=True)
>>> print(res)
LipinskiRule(MW=228.29, logP=5.15, nHD=0, nHA=0, Disposed='Accept', nViolated=1)
```

If you are  not familiar with some rules, you could use <code>help</code>function

```python
>>> help(CheckRule.CheckEganRule)
Help on function CheckEganRule in module scopy.Druglikeness.CheckRule:

CheckEganRule(mol, detail)
    #################################################################
    Bad or Good oral biovailability rule
    
    -Ref.:
        Egan, William J., Kenneth M. Merz, and John J. Baldwin. 
        Journal of medicinal chemistry 43.21 (2000): 3867-3877.
        
    -Rule details:
        0 <= tPSA <= 132
        -1 <= logP <=6
    #################################################################
```

you could also customize your rules

```python
>>> res = CheckRule.Check_CustomizeRule(mol,prop_kws={'MW':(None,500),
                                                  'nRing':(1,5),
                                                  'nHB':(1,5)},
                                        	detail=True,closed_interval=False)
>>> print(res)
CustomizeRule(MW=228.29, nRing=4, nHB=0, nViolate=1, VioProp=['nHB'])
```

### Filter compounds through SMARTS

The filters consist of a series of molecular query strings written using the SMARTS coding language described by [Daylight](https://www.daylight.com/). 

```python
"""
---
Up to now(2019.07.02), we have collected followed endpoints(the number of SMARTS):
	Acute_Aquatic_Toxicity(99)
	AlphaScreen_FHs(6)
	AlphaScreen_GST_FHs(34)
	AlphaScreen_HIS_FHs(19)
	Biodegradable(9)
	BMS(180)
	Chelating(55)
	Developmental_Mitochondrial(12)
	Genotoxic_Carcinogenicity_Mutagenicity(117)
	Idiosyncratic(35)
	LD50_oral(20)
	Luciferase_Inhibitory(3)
	NonBiodegradable(19)
	NonGenotoxic_Carcinogenicity(23)
	NTD(105)
	Pains(480)
	Potential_Electrophilic(119)
	Promiscuity(177)
	Reactive_Unstable_Toxic(335)
	Skin_Sensitization(155)
	SureChEMBL(165)
---	
Total: 23 endpints with 2167 SMARTS
"""
```

Besides, we have collected **450 SMARTS** about **E**xtended **F**unctional **G**roups(EFG), an efficient set
for chemical characterization and structure-activity relationship studies of chemical compounds.

```python
>>> from scopy.StructureAlert import FliterWithSmarts
>>> mol = Chem.MolFromSmiles('C1=CC=CC2C=CC3C4C=CC=CC=4C=CC=3C1=2')
>>> res = FliterWithSmarts.Check_SureChEMBL(mol)
print(res)
CheckRes(Disposed='Rejected', Endpoint='SureChEMBL')
```

Similarly, you could pass <code>True</code> to parameter show to get more information.

```python
>>> res = FliterWithSmarts.Check_SureChEMBL(mol,detail=True)
CheckRes(Disposed='Rejected', MatchedAtoms=[((3, 2, 1, 0, 17, 16, 15, 14, 13, 8, 7, 6, 5, 4), (12, 11, 10, 9, 8, 7, 6, 5, 4, 17, 16, 15, 14, 13))], MatchedNames=['Polynuclear_Aromatic_2'], Endpoint='SureChEMBL')
```

### Visualization

