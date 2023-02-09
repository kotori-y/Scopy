<!--
 * @Description: 
 * @Author: Kotori Y
 * @Date: 2019-09-15 17:04:40
 * @LastEditors: Kotori Y
 * @LastEditTime: 2021-05-15 21:40:28
 * @FilePath: \scopy\README.md
 * @AuthorMail: kotori@cbdd.me
-->
# Scopy: An integrated negative design Python library for desirable HTS/VS database design

[![Travis (.com)](https://img.shields.io/travis/com/kotori-y/scopy)](https://travis-ci.com/kotori-y/Scopy) [![codecov](https://codecov.io/gh/kotori-y/Scopy/branch/master/graph/badge.svg)](https://codecov.io/gh/kotori-y/Scopy) [![GitHub last commit](https://img.shields.io/github/last-commit/kotori-y/scopy)](https://github.com/kotori-y/Scopy/commits/master) [![Conda](https://img.shields.io/badge/Install%20with-conda-green)](https://conda.anaconda.org/kotori_y) [![PyPI](https://img.shields.io/badge/Install%20with-pypi-informational)](https://pypi.org/project/scopy/) [![MIT License](https://img.shields.io/badge/license-MIT-black)](https://anaconda.org/kotori_y/scopy) [![Blog](https://img.shields.io/badge/blog-iamkotori-pink)](https://blog.kotori.wiki/) [![Kouhai](https://img.shields.io/badge/contributor-Ziyi-%23B3D0BE)](https://github.com/Yangziyi1997) [![DOI](https://img.shields.io/badge/doi-Briefings%20in%20Bioinformatics-informational)](https://doi.org/10.1093/bib/bbaa194)

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

## Cite us

Yang ZY, Yang ZJ, Lu AP, Hou TJ, Cao DS. Scopy: an integrated negative design python library for desirable HTS/VS database design [published online ahead of print, 2020 Sep 7]. *Brief Bioinform*. 2020;bbaa194. doi:10.1093/bib/bbaa194

```
@article{10.1093/bib/bbaa194,
    author = {Yang, Zi-Yi and Yang, Zhi-Jiang and Lu, Ai-Ping and Hou, Ting-Jun and Cao, Dong-Sheng},
    title = "{Scopy: an integrated negative design python library for desirable HTS/VS database design}",
    journal = {Briefings in Bioinformatics},
    year = {2020},
    month = {09},
    abstract = "{High-throughput screening (HTS) and virtual screening (VS) have been widely used to identify potential hits from large chemical libraries. However, the frequent occurrence of ‘noisy compounds’ in the screened libraries, such as compounds with poor drug-likeness, poor selectivity or potential toxicity, has greatly weakened the enrichment capability of HTS and VS campaigns. Therefore, the development of comprehensive and credible tools to detect noisy compounds from chemical libraries is urgently needed in early stages of drug discovery.In this study, we developed a freely available integrated python library for negative design, called Scopy, which supports the functions of data preparation, calculation of descriptors, scaffolds and screening filters, and data visualization. The current version of Scopy can calculate 39 basic molecular properties, 3 comprehensive molecular evaluation scores, 2 types of molecular scaffolds, 6 types of substructure descriptors and 2 types of fingerprints. A number of important screening rules are also provided by Scopy, including 15 drug-likeness rules (13 drug-likeness rules and 2 building block rules), 8 frequent hitter rules (four assay interference substructure filters and four promiscuous compound substructure filters), and 11 toxicophore filters (five human-related toxicity substructure filters, three environment-related toxicity substructure filters and three comprehensive toxicity substructure filters). Moreover, this library supports four different visualization functions to help users to gain a better understanding of the screened data, including basic feature radar chart, feature-feature-related scatter diagram, functional group marker gram and cloud gram.Scopy provides a comprehensive Python package to filter out compounds with undesirable properties or substructures, which will benefit the design of high-quality chemical libraries for drug design and discovery. It is freely available at https://github.com/kotori-y/Scopy.}",
    issn = {1477-4054},
    doi = {10.1093/bib/bbaa194},
    url = {https://doi.org/10.1093/bib/bbaa194},
    note = {bbaa194},
    eprint = {https://academic.oup.com/bib/advance-article-pdf/doi/10.1093/bib/bbaa194/33719387/bbaa194.pdf},
}
```

  

## Acknowledgement

*Thanks to my colleague, [Ziyi](https://github.com/Yangziyi1997), for assisting me to complete the writing of document and article.*
