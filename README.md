# Scopy: a compounds filter for HTS and VS

[![Build Status](https://travis-ci.com/kotori-y/Scopy.svg?branch=master)](https://travis-ci.com/kotori-y/Scopy) [![Read the Docs](https://img.shields.io/readthedocs/scopy)](https://scopy.readthedocs.io/en/latest/) [![GitHub last commit](https://img.shields.io/github/last-commit/kotori-y/scopy)](https://github.com/kotori-y/Scopy/commits/master) [![Blog](https://img.shields.io/badge/blog-iamkotori-informational)](https://blog.iamkotori.com/) [![License](https://img.shields.io/github/license/kotori-y/scopy)](https://opensource.org/licenses/MIT) [![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu)

<div align=center>
    <img src='Scopy.png'>
</div>

## Overview

**Scopy(<font color='red'>S</font>creening <font color='red'>CO</font>mpounds)**, based on RDKit, is an integrated negative design python library designed for screening out undesiable compounds in the early drug discovery. Scopy includes six modules, covering **data preparation**, **screening filters**, the **calculation of scaffolds and descriptors**, and the **visualization analysis**.

## Installation

### Install RDKit

```
>>> conda install -c conda-forge rdkit
```

### Install Scopy

Scopy has been successfully tested on Linux and Windows systems under python3 enviroment.

```
>>> git clone git@github.com:kotori-y/Scopy.git && cd scopy
>>> [sudo] python setup.py install
```

## Documentation

(1)The online version of the documentation is available here: https://scopy.readthedocs.io/en/latest/<br>(2)Quick start examples: https://scopy.readthedocs.io/en/latest/user_guide.html<br>(3)Application examples(pipelines): https://scopy.readthedocs.io/en/latest/application.html

## Contact

If you have questions or suggestions, please contact: kotori@cbdd.me,and oriental-cds@163.com.

Please see the file LICENSE for details about the "MIT" license which covers this software and its associated data and documents.
