# Scopy: An integrated negative design Python library for desirable HTS/VS database design

[![Travis (.com)](https://img.shields.io/travis/com/kotori-y/scopy)](https://travis-ci.com/kotori-y/Scopy) [![codecov](https://codecov.io/gh/kotori-y/Scopy/branch/master/graph/badge.svg)](https://codecov.io/gh/kotori-y/Scopy) [![GitHub last commit](https://img.shields.io/github/last-commit/kotori-y/scopy)](https://github.com/kotori-y/Scopy/commits/master) [![Conda](https://img.shields.io/badge/Install%20with-conda-green)](https://conda.anaconda.org/kotori_y) [![PyPI](https://img.shields.io/badge/Install%20with-pypi-informational)](https://pypi.org/project/scopy/) [![MIT License](https://img.shields.io/badge/license-MIT-black)](https://anaconda.org/kotori_y/scopy) [![Blog](https://img.shields.io/badge/blog-iamkotori-pink)](https://blog.iamkotori.com/) [![Kouhai](https://img.shields.io/badge/contributor-Ziyi-%23B3D0BE)](https://github.com/Yangziyi1997)

<div align=center>
    <img src='Scopy.png'>
</div>

## Overview

**Scopy** (**S**crenning **CO**mpounds in **PY**thon), an integrated negative design python library designed for screening out undesirable compounds in the early drug discovery. Scopy includes six modules, **covering data preparation**, **screening filters**, the **calculation of scaffolds** and **descriptors**, and the **visualization analysis**. 

## Installation

### Install RDKit

```
>>> conda install -c conda-forge rdkit
```

### Install Scopy

Scopy has been successfully tested on Linux, OSX and Windows systems under Python3 enviroment.

#### Source

```
>>> git clone git@github.com:kotori-y/Scopy.git && cd scopy
>>> [sudo] python setup.py install
```

#### Conda [![Conda](https://img.shields.io/conda/v/kotori_y/scopy?color=green&label=conda&style=flat-square)](https://anaconda.org/kotori_y/scopy)

```
>>> conda install -c kotori_y scopy
```

#### PyPI [![PyPI](https://img.shields.io/pypi/v/scopy?style=flat-square)](https://pypi.org/project/scopy/)

```
>>> pip install scopy
```

## Documentation

(1)The online version of the documentation is available here: https://scopy.iamkotori.com/<br>(2)Quick start examples: https://scopy.iamkotori.com/user_guide.html<br>(3)Application examples(pipelines): https://scopy.iamkotori.com/application.html

## Contact

If you have questions or suggestions, please contact: kotori@cbdd.me,and oriental-cds@163.com.<br>Please see the file LICENSE for details about the "MIT" license which covers this software and its associated data and documents.

## Acknowledgement

*Thanks to my colleague, [Ziyi](https://github.com/Yangziyi1997), for assisting me to complete the writing of document and article.*