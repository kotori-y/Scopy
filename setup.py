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

# print（__doc__）

package_data= {'scopy':['structure_alert/*','druglikeness/*','fingerprint/*',
                        'visualize/*','pretreat/*','data/SMARTS/*','data/PATT/*',
                        'data/ACID/*','data/*','data/Crippen/*','data/EFG/*','data/Demo/*',
                        'data/MC/mcloud/*','data/MC/mcloud/ertl/mcloud/*',
                        'data/MOL/*','data/PubChem/*',]}

#package_data= {'scopy':['structure_alert/*','druglikeness/*','test/*','data/SMARTS/*','data/PATT/*','data/ACID/*','data/*','data/Crippen/*','fingerprint/*']}


setup(name="cbdd-scopy",  
      version="1.1.2",
      license="MIT",
      description="A filter tool for HTS and VS",
      long_description="Scopy (Screening COmpounds in PYthon), based on RDKit, is an integrated negative design python library designed for screening out undesiable compounds in the early drug discovery.",
      author="Zhi-Jiang Yang (Kotori), Dong-Sheng Cao",
      author_email="yzjkid9@gmail.com",
      maintainer="Zhi-Jiang Yang (Kotori)",
      maintainer_email="kotori@cbdd.me",
      url="https://github.com/kotori-y/Scopy",
      package_data=package_data,
      # include_package_data=True,
      package_dir={'scopy':'scopy'},
      py_modules = ['scopy.ScoConfig'],
      packages=['scopy']
      )
