# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 09:46:11 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

"""

__doc__="""
        Words
        """

from scopy.Druglikeness import rulefilter
from rdkit import Chem
import matplotlib.pyplot as plt
import numpy as np

Rule_dict = {
        'BeyondRo5':{'MW':(0,1000),'logP':(-2,10),'nHD':(0,6),'nHA':(0,15),'TPSA':(0,250),'nRot':(0,20)},
        'Egan':{'tPSA':(0,132),'logP':(-1,6)},
        'Veber':{'nRot':(0,10),'tPSA':(0,140),'nHB':(0,12)},
        'Lipinski':{'MW':(0,500),'logP':(0,5),'nHD':(0,5),'nHA':(0,10)},
        'Xu':{'nHD':(0,5),'nHA':(0,10),'nRot':(3,35),'nRing':(1,7),'nHev':(10,50)}
        }



def _Radar(CheckResult,Rule):
    res = [x for x in CheckResult[:-2]]
    length = len(res)
    
    rule_ceil = np.array([x[1] for x in Rule_dict[Rule].values()])
    rule_floor = np.array([x[0] for x in Rule_dict[Rule].values()])
    
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
    ax.set_xticklabels([x for x in Rule_dict[Rule].keys()]) 
    ax.set_yticks([])
    
    ax.spines['polar'].set_visible(False)
    ax.grid(axis='y')
    ax.set_ylim([0,1])
    for i in [0,0.2,0.4,0.6,0.8,1.0]:
        ax.plot(X_ticks,[i]*(length+1),'-', color='black',lw=0.5)
    ax.set_theta_zero_location('N')
#    plt.show()
    
    
def VisualizeBeyondRo5(mol=None,bRo5Rule=None):
    if not mol:
        _Radar(bRo5Rule,'BeyondRo5')
    else:
        bRo5Rule = rulefilter.CheckBeyondRo5(mol,detail=True)
        _Radar(bRo5Rule,'BeyondRo5')
    

def VisualizeLipinski(mol=None,LipinskiRule=None):
    if not mol:
        _Radar(LipinskiRule,'Lipinski')
    else:
        LipinskiRule = rulefilter.CheckLipinskiRule(mol,detail=True)
        _Radar(LipinskiRule,'Lipinski')


def VisualizeXu(mol=None,XuRule=None):
    if not mol:
        _Radar(XuRule,'Xu')
    else:
        XuRule = rulefilter.CheckXuRule(mol,detail=True)
        _Radar(XuRule,'Xu')

def VisualizeEgan(EganRule=None,mol=None):
    if mol:
        EganRule = rulefilter.CheckEganRule(mol)
    else:
        pass
    y,x = [x for x in EganRule[:-2]]
    f,ax = plt.subplots()
    
    y_ceil = max(y,132)
    ax.fill([-1,6,2.5,-1],[0,0,132,0],'#D7FCCB',
            zorder=0)
    ax.scatter(x,y)
#    ax.set_xlim(0)
    ax.set_ylim([0,1.1*y_ceil])
    ax.tick_params(right=True,top=True,width=1,direction='in')
    ax.set_xlabel('LogP')
    ax.set_ylabel('tPSA')
    plt.show()



def PfizerPositioning(PfizerRule):
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
    
    


if '__main__' == __name__:
    mol = Chem.MolFromSmiles('Fc1ccc(CC2=NNC(=O)c3ccccc23)cc1C(=O)N4CCc5cccc6C(=O)NCC4c56')
#    bRo5Rule = rulefilter.CheckBeyondRo5(mol,detail=True)
##    VisualizeBeyondRo5(bRo5Rule)
##    PfizerRule = CheckRule.CheckPfizerRule(mol,detail=True)
##    PfizerPositioning(PfizerRule)
#    
#    VeberRule = rulefilter.CheckVeberRule(mol,detail=True)
    Egan = rulefilter.CheckEganRule(mol,detail=True)
#    VisualizeLipinski(mol=mol)
    VisualizeEgan(Egan)