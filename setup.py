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



<<<<<<< HEAD
package_data= {'scopy':['structure_alert/*','druglikeness/*','test/*','data/SMARTS/*','data/PATT/*','data/ACID/*','data/*','data/Crippen/*','data/EFG/*','fingerprint/*','visualize/*','pretreat/*']}
=======
package_data= {'scopy':['structure_alert/*','druglikeness/*','test/*','data/SMARTS/*','data/PATT/*','data/ACID/*','data/*','data/Crippen/*','fingerprint/*']}
>>>>>>> d9ae338fbdfa8b3e677463b33d5e8237a5528729

setup(name="scopy",  
      version="1.1", 
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
