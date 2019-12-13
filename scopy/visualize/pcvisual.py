# -*- coding: utf-8 -*-

#Created on Mon Sep 16 10:25:57 2019
#
#@Author: Zhi-Jiang Yang, Dong-Sheng Cao
#@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
#@Homepage: http://www.scbdd.com
#@Mail: yzjkid9@gmail.com; oriental-cds@163.com
#@Blog: https://blog.moyule.me


import matplotlib.pyplot as plt
import numpy as np
try:
    from ..druglikeness import rulesfilter
    from ..druglikeness import molproperty
except:
    import sys
    sys.path.append('..')
    from druglikeness import rulesfilter
    from druglikeness import molproperty



class Propvisual(object):
    """
    """
    def __init__(self, mol):
        self.mol = mol
        self.Rule_dict = {
        'BeyondRo5':{'MW':(0,1000),'logP':(-2,10),'nHD':(0,6),'nHA':(0,15),'TPSA':(0,250),'nRot':(0,20)},
        'Egan':{'tPSA':(0,132),'logP':(-1,6)},
        'Veber':{'nRot':(0,10),'tPSA':(0,140),'nHB':(0,12)},
        'Lipinski':{'MW':(0,500),'logP':(0,5),'nHD':(0,5),'nHA':(0,10)},
        'Xu':{'nHD':(0,5),'nHA':(0,10),'nRot':(3,35),'nRing':(1,7),'nHev':(10,50)}
        }
 
    def _Radar(self, CheckResult, Rule):
        res = [x for x in CheckResult[:-2]]
        length = len(res)
        
        rule_ceil = np.array([x[1] for x in self.Rule_dict[Rule].values()])
        rule_floor = np.array([x[0] for x in self.Rule_dict[Rule].values()])
        
        bench_floor = np.vstack((res,rule_floor)).min(axis=0)
        bench_ceil = np.vstack((res,rule_ceil)).max(axis=0)*1.2
        
        res = (res-bench_floor)/(bench_ceil-bench_floor)
        floor = (rule_floor-bench_floor)/(bench_ceil-bench_floor)
        ceil = (rule_ceil-bench_floor)/(bench_ceil-bench_floor)
        
        theta = np.linspace(0, 360, length, endpoint=False)
        X_ticks = np.radians(theta)#angle to radian
        X_ticks = np.append(X_ticks,X_ticks[0])
        
        Y = np.vstack((res,floor,ceil))
        Y = np.hstack((Y, Y[:,0].reshape(3,1)))
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.plot(X_ticks, Y[0])
        ax.plot(X_ticks, Y[1],color='#FF7B9A')
        ax.fill(X_ticks, Y[1], alpha=0.25,color='#FF7B9A')
        ax.plot(X_ticks, Y[2],color='#66B9EF')
        ax.fill(X_ticks, Y[2], alpha=0.25,color='#66B9EF')
        ax.set_xticks(X_ticks)
        ax.set_xticklabels([x for x in self.Rule_dict[Rule].keys()]) 
        ax.set_yticks([])
        
        ax.spines['polar'].set_visible(False)
        ax.grid(axis='y')
        ax.set_ylim([0,1])
        for i in [0,0.2,0.4,0.6,0.8,1.0]:
            ax.plot(X_ticks,[i]*(length+1),'-', color='black',lw=0.5)
        ax.set_theta_zero_location('N')
        plt.show()
        
    def VisualizeBeyondRo5(self):
        bRo5Rule = rulesfilter.CheckBeyondRo5(self.mol, detail=True)
        self._Radar(bRo5Rule,'BeyondRo5')
        
    def VisualizeLipinski(self):
        LipinskiRule = rulesfilter.CheckLipinskiRule(self.mol,detail=True)
        self._Radar(LipinskiRule,'Lipinski')
            



#def propdist(mols=None, pair=False, pos_mols=None, neg_mols=None):
#    def getprops(mols):
#        logPs = list(map(lambda mol: molproperty.CalculateLogP(mol), mols))
#        TPSAs = list(map(lambda mol: molproperty.CalculateTPSA(mol), mols))
#        mws = list(map(lambda mol: molproperty.CalculateMolWeight(mol), mols))
#        nrot = list(map(lambda mol: molproperty.CalculateNumRotatableBonds(mol), mols))
#        nHD = list(map(lambda mol: molproperty.CalculateNumHDonors(mol), mols))
#        nHA = list(map(lambda mol: molproperty.CalculateNumHAcceptors(mol), mols))
#        props = [logPs,TPSAs,mws,nrot,nHD,nHA]
#        return props
#    
#    labels = ['logP', 'TPSA', 'MW', 'nRot', 'nHD', 'nHA']
#    if not pair:
#        props = getprops(mols)
#        colors = ['#814c94', '#4884ff', '#9a54f7', '#eb50ff', '#fb70a2', '#f7104a'
#                  ]
#        
#        
#        f,axes = plt.subplots(2,3,figsize=(24,10))
#        for ax,prop,color,label in zip(axes.flatten(), props, colors, labels):
#            ax.hist(prop,ec='black',fc=color,alpha=0.7)
#            ax.spines['right'].set_color('None')
#            ax.spines['top'].set_color('None')
#            ax.spines['bottom'].set_linewidth(2)
#            ax.spines['left'].set_linewidth(2)
#            ax.set_xlabel(label, fontdict={'family':'arial','size':18,'color':color})
#        plt.show()
#    else:
#        pos_colors = ['#814c94', '#4884ff', '#9a54f7', '#eb50ff', '#fb70a2', '#f7104a'
#                  ]
#        neg_colors = ['#5f944c', '#ffc348', '#b1f754', '#64ff50', '#70fbc9', '#10f7bd']
#        pos_props = getprops(pos_mols)
#        neg_props = getprops(neg_mols)
#        
#        f,axes = plt.subplots(2,3,figsize=(24,10))
#        for ax,pos_prop,pos_color,neg_prop,neg_color,label in zip(axes.flatten(), pos_props, pos_colors, neg_props, neg_colors, labels):
#            ax.hist(pos_prop,ec='black',fc=pos_color,alpha=0.6,label='Pos')
#            ax.hist(neg_prop,ec='black',fc=neg_color,alpha=0.6,label='Neg')
#            ax.spines['right'].set_color('None')
#            ax.spines['top'].set_color('None')
#            ax.spines['bottom'].set_linewidth(2)
#            ax.spines['left'].set_linewidth(2)
#            ax.set_xlabel(label, fontdict={'family':'arial','size':18})
#            ax.legend()
#        plt.show()
        

def propdist(mols=None):
    from itertools import product
    def getprops(mols):
        logPs = list(map(lambda mol: molproperty.CalculateLogP(mol), mols))
        TPSAs = list(map(lambda mol: molproperty.CalculateTPSA(mol), mols))
        mws = list(map(lambda mol: molproperty.CalculateMolWeight(mol), mols))
        nrot = list(map(lambda mol: molproperty.CalculateNumRotatableBonds(mol), mols))
        nHD = list(map(lambda mol: molproperty.CalculateNumHDonors(mol), mols))
        nHA = list(map(lambda mol: molproperty.CalculateNumHAcceptors(mol), mols))
        props = dict(zip(['logP','TPSA','MW','nRot','nHD','nHA'],
                         [logPs,TPSAs,mws,nrot,nHD,nHA]))
        return props
    

    f,axes = plt.subplots(6,6,figsize=(25,25))
    labels = ['logP', 'TPSA', 'MW', 'nRot', 'nHD', 'nHA']
    props = getprops(mols)
    
    for pair in product(range(6), repeat=2):
        col,row = pair[0],pair[1]
        x_label = labels[row]
        y_label = labels[col]
        axes[col][row].spines['right'].set_color('None')
        axes[col][row].spines['top'].set_color('None')
        axes[col][row].spines['bottom'].set_linewidth(2)
        axes[col][row].spines['left'].set_linewidth(2)
        axes[col][row].tick_params(length=8)
        if col == row:
            axes[col][row].hist(props[x_label], ec='black',fc='#66cdaa')
        else:
            axes[col][row].scatter(props[x_label], props[y_label], edgecolors='black', alpha=0.8, color='#66cdaa')
        if col == 5:
            axes[col][row].set_xlabel(x_label, fontdict={'family':'arial','size':16})
        else:
            axes[col][row].set_xticklabels([])
        if row == 0:
            axes[col][row].set_ylabel(y_label, fontdict={'family':'arial','size':16})
        else:
            axes[col][row].set_yticklabels([])
    plt.show()


 
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
    








#def piegroup(mols):
#    fgs = map(lambda mol: molproperty.GetIFG(mol), mols)
#    fgs = [i.type for fg in fgs for i in fg]
#    return fgs
    
    
    
if '__main__' == __name__:
    from rdkit.Chem import AllChem as Chem
    mol = Chem.MolFromSmiles('Fc1ccc(CC2=NNC(=O)c3ccccc23)cc1C(=O)N4CCc5cccc6C(=O)NCC4c56')
#    neg_mols = Chem.SDMolSupplier(r"C:\Users\0720\Documents\Tencent Files\1223821976\FileRecv\ace_dud_decoys.sdf")
#    pos_mols = Chem.SDMolSupplier(r"C:\Users\0720\Documents\Tencent Files\1223821976\FileRecv\actives_hsp90.sdf")
#    propdist(pair=True,pos_mols=pos_mols,neg_mols=neg_mols)
    v = Propvisual(mol)
    v.VisualizeBeyondRo5()
#    
    

    
        