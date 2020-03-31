# Scopy: a compounds filter for HTS and VS

[![Travis (.com)](https://img.shields.io/travis/com/kotori-y/scopy?style=flat-square)](https://travis-ci.com/kotori-y/Scopy) [![Read the Docs](https://img.shields.io/readthedocs/scopy?style=flat-square)](https://scopy.readthedocs.io/en/latest/) [![GitHub last commit](https://img.shields.io/github/last-commit/kotori-y/scopy?style=flat-square)](https://github.com/kotori-y/Scopy/commits/master) [![Conda](https://anaconda.org/kotori_y/scopy/badges/installer/conda.svg)](https://conda.anaconda.org/kotori_y) [![PyPI](https://img.shields.io/badge/Install%20with-pypi-informational?style=flat-square)](https://pypi.org/project/cbdd-scopy/) [![MIT License](https://anaconda.org/kotori_y/scopy/badges/license.svg)](https://anaconda.org/kotori_y/scopy) [![Blog](https://img.shields.io/badge/blog-iamkotori-pink?style=flat-square)](https://blog.iamkotori.com/) [![996.icu](https://img.shields.io/badge/link-996.icu-%23FF4D5B.svg?style=flat-square)](https://996.icu/#/en_US)

<div align=center>
    <img src='Scopy.png'>
</div>

## Overview

**Scopy(<font color='red'>S</font>creening <font color='red'>CO</font>mpounds in <font color='red'>PY</font>thon)**, based on RDKit, is an integrated negative design python library designed for screening out undesiable compounds in the early drug discovery. Scopy includes six modules, covering **data preparation**, **screening filters**, the **calculation of scaffolds and descriptors**, and the **visualization analysis**.

## Installation

### Install RDKit

```
>>> conda install -c conda-forge rdkit
```

### Install Scopy

Scopy has been successfully tested on Linux and Windows systems under python3 enviroment.

#### Source

```
>>> git clone git@github.com:kotori-y/Scopy.git && cd scopy
>>> [sudo] python setup.py install
```

#### Conda

```
>>> conda install kotori_y scopy
```

#### PyPI

```
>>> pip install cbdd-scopy
```

## Documentation

(1)The online version of the documentation is available here: https://scopy.readthedocs.io/en/latest/<br>(2)Quick start examples: https://scopy.readthedocs.io/en/latest/user_guide.html<br>(3)Application examples(pipelines): https://scopy.readthedocs.io/en/latest/application.html

## Contact

If you have questions or suggestions, please contact: kotori@cbdd.me,and oriental-cds@163.com.<br>Please see the file LICENSE for details about the "MIT" license which covers this software and its associated data and documents.

