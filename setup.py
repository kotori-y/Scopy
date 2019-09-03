# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:21:07 2019

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

I love my senpai forerver!:P
"""

from __future__ import absolute_import
from distutils.core import setup



package_data= {'scopy':['StructureAlert/*','Druglikeness/*','test/*','data/SMARTS/*','data/PATT/*','data/ACID/*','data/*','data/Crippen/*','fingerprint/*']}

setup(name="scopy",  
      version="1.0", 
      description="",
      long_description="",
      author="Orient&Kotori_Y",
      author_email="yzjkid9@gmail.com",
      url="blog.moyule.me",
      package_data=package_data,
      package_dir={'scopy':'scopy'},
      py_modules = ['scopy.ScoConfig'],
      packages=['scopy']
      )  
