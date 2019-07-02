# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:21:07 2019

You are not expected to understand my codes!

@Author: Kotori_Y
@Blog: blog.moyule.me
@Weibo: Kotori-Y
@Mail: yzjkid9@gmial.com

I love Megumi forerver!
"""

from __future__ import absolute_import
from distutils.core import setup

package_data= {'scopy':['StructureAlert/*','Druglikeness/*','test/*','data/SMARTS/*','data/*']}

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
