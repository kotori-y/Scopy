# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 15:08:41 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

♥I love Princess Zelda forever♥
"""

import shutil
from subprocess import run

def ShowMcloud(file, number=150, skip=0, savedir=None, hidden=False):
    """
    """
    i = '-i' if savedir else ''
    nogui = '-nogui' if hidden else ''
    
    command = 'cd /d .\data\mcloud && \
    java -cp ".;depictjni.jar;depict.jar" ertl/mcloud/MCloud\
    -f {} -n {} -skip {} {} {}'.format(file, number, skip, i, nogui)
          
    run(command,shell=True)
    
    try:
        shutil.move('.\data\mcloud\mcloud.png',savedir)
    except FileNotFoundError:
        pass
    
if '__main__'==__name__:
    ShowMcloud(r"scaffolds.smi",savedir=r"mcloud.png")