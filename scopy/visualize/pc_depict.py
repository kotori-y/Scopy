# -*- coding: utf-8 -*-

#Created on Fri Dec  6 16:10:33 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me



from itertools import product
import matplotlib.pyplot as plt
import numpy as np

try:
    from ..druglikeness import rulesfilter
    from ..druglikeness import molproperty
    from ..druglikeness import druglikeness
except:
    import sys
    sys.path.append('..')
    from druglikeness import rulesfilter
    from druglikeness import molproperty
    from druglikeness import druglikeness
    


def prop_matrix(mols,n_jobs=1,items=['logP','TPSA','MW','nRot','nHD','nHA']):
    """The proprty matrix can intuitively show the compounds' distribution 
    in Two-Dimension space, and diagonal of the matrix is the displot of property
    
    :param mols: Molecules
    :type mols: iterable object
    :param n_jobs: The number of CPUs to use to do the computation, defaults to 1
    :type n_jobs: int, optional
    :param items: The abbreviation name of property to be visualized, defaults to ['logP','TPSA','MW','nRot','nHD','nHA']
    :type items: list, optional
    :return: The proerty matrix
    :rtype: matplotlib.figure.Figure
        
    """
    props = druglikeness.PC_properties(mols,n_jobs)
    props = props.GetProperties(items=items)
    length = len(items)
    
    fig,axes = plt.subplots(length,length,figsize=(4*length,4*length))
    
    for pair in product(range(len(items)), repeat=2):
        col,row = pair[0],pair[1]
        x_label = items[row]
        y_label = items[col]
        axes[col][row].spines['right'].set_color('None')
        axes[col][row].spines['top'].set_color('None')
        axes[col][row].spines['bottom'].set_linewidth(2)
        axes[col][row].spines['left'].set_linewidth(2)
        axes[col][row].tick_params(length=2)
        axes[col][row].tick_params(direction='in', which='both')
        if col == row:
            axes[col][row].hist(props[x_label], ec='black',fc='#DA3B2B')
        else:
            axes[col][row].scatter(props[x_label], props[y_label], 
                edgecolors='black', alpha=0.8, s=23, color='#DA3B2B')
        if col == length-1:
            axes[col][row].set_xlabel(x_label, fontdict={'family':'arial','size':16})
        else:
            axes[col][row].set_xticklabels([])
        if row == 0:
            axes[col][row].set_ylabel(y_label, fontdict={'family':'arial','size':16})
        else:
            axes[col][row].set_yticklabels([])
    
    plt.show()
    return fig
  
          
def rule_radar(mol, 
               prop_kws={'MW':(100,600),'logP':(-3,6),
                         'nHA':(0,12),'nHD':(0,7),'TPSA':(0,180),
                         'nRot':(0,11),'nRing':(0,6),'MaxRing':(0,18),
                         'nC':(3,35),'nHet':(1,15),'HetRatio':(0.1,1.1),
                         'fChar':(-4,4),'nRig':(0,30)}
               ):
    """ A radar plot positionning compound's values within the selected filter ranges (pale blue and red). 
    By default, the `drug-like soft`_ filter ranges are visualized.
    
    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: A radar plot positionning compound's values
    :rtype: matplotlib.figure.Figure
    
    .. _drug-like soft:
        http://fafdrugs4.mti.univ-paris-diderot.fr/filters.html

    """
    def _disposedNone(dtype,num_a,num_b):
        if dtype == 'MIN':
            NUM = min([num_a,num_b])
            return NUM - 1.2*abs(NUM)
        else:
            NUM = max([num_a,num_b])
            return NUM + 1.2*abs(NUM)
        
    items = list(prop_kws.keys())
    props = molproperty.GetProperties(mol,items=items)
    num_prop = len(props)
    
    for item in items:
        assert (prop_kws[item][0] is not None or 
                prop_kws[item][1] is not None
                ), "You need to enter at least an upper or lower limit"
        
    rule_ceil = np.array([prop_kws[item][1] 
    if prop_kws[item][1] is not None else _disposedNone('MAX',prop_kws[item][0],props[item]) 
    for item in items])
    
    rule_floor = np.array([prop_kws[item][0]
    if prop_kws[item][0] is not None else _disposedNone('MIN',prop_kws[item][1],props[item]) 
    for item in items])
    
    props = np.array(list(props.values()))

    bench_floor = np.vstack((props, rule_floor)).min(axis=0)
    bench_floor -= 0.2*bench_floor
    bench_ceil = np.vstack((props, rule_ceil)).max(axis=0)*1.2
    
    #max-min standarize
    props = (props-bench_floor)/(bench_ceil-bench_floor)
    floor = (rule_floor-bench_floor)/(bench_ceil-bench_floor)
    ceil = (rule_ceil-bench_floor)/(bench_ceil-bench_floor)                 
    
    theta = np.linspace(0, 360, num_prop, endpoint=False)
    X_ticks = np.radians(theta)#angle to radian
    X_ticks = np.append(X_ticks,X_ticks[0])
    Y = np.vstack((props,floor,ceil))
    Y = np.hstack((Y, Y[:,0].reshape(3,1)))
    
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.plot(X_ticks, Y[0])
    ax.plot(X_ticks, Y[1],color='#FF7B9A')
    ax.fill(X_ticks, Y[1], alpha=0.25,color='#FF7B9A')
    ax.plot(X_ticks, Y[2],color='#EDB035')
    ax.fill(X_ticks, Y[2], alpha=0.20,color='#EDB035')
    ax.set_xticks(X_ticks)
    ax.set_xticklabels(items) 
    ax.set_yticks([])
    
    ax.spines['polar'].set_visible(False)
    ax.grid(axis='y')
    ax.set_ylim([0,1])
    for i in [0,0.2,0.4,0.6,0.8,1.0]:
        ax.plot(X_ticks,[i]*(num_prop+1),'-', color='black',lw=0.5)
    ax.set_theta_zero_location('N')
    
    plt.show()           
    return fig
    
    
def oral_absorption_radar(mol):
    """
    """
    prop_kws={'logP':(-2,5),'MW':(150,500),'TPSA':(20,150),
              'nRot':(0,10),'nHA':(0,10),'nHD':(0,5)}
    
    fig = rule_radar(mol, prop_kws)

    return fig
    
    
def PfizerPositioning(mol=None, PfizerRule=None):
    if mol:
        PfizerRule = rulesfilter.CheckPfizerRule(mol, detail=True)
    else:
        pass
    f,ax = plt.subplots()
    res = [x for x in PfizerRule[:-2]]
    ax.scatter(*res,s=30,color='black')
    
    ax.fill([-10,3,3,-10,-10],[0,0,75,75,0],'#D7FCCB',
            [3,10,10,3,3],[0,0,75,75,0],'#FF7B9A',
            [3,10,10,3,3],[75,75,200,200,75],'#D7FCCB',
            [-10,3,3,-10,-10],[75,75,200,200,75],'#93FFA3',zorder=0)    
    ax.set_xticks([-10,3,10])
    ax.set_yticks([0,75,200])
    ax.set_xlim([-10,10])
    ax.set_ylim([0,200])
    ax.set_xlabel('LogP')
    ax.set_ylabel('tPSA')
    plt.show()
    
if '__main__'==__name__:
    from rdkit import Chem
    
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
            'CCCCC(=O)[H]', #Biodegradable
            'C1=CN=C(C(=O)O)C=C1', #Chelating
            'C(OC)1=CC=C2OCC3OC4C=C(OC)C=CC=4C(=O)C3C2=C1',
            'C1=C2N=CC=NC2=C2N=CNC2=C1', #Genotoxic_Carcinogenicity_Mutagenicity
            'N(CC)(CCCCC)C(=S)N', #Idiosyncratic
            ]
    mols = (Chem.MolFromSmiles(smi) for smi in smis)
    fig = prop_matrix(mols)
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    