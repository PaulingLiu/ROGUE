# ROGUE (Ratio of Global Unshifted Entropy)

## Contents

- [Overview](#overview)
- [Installation Guide](#installation-guide)
- [Tutorial](#tutorial)
- [License](./LICENSE)
- [Citation](#citation)
- [Contact](#Contact)

## Overview
Often, it is not even clear whether a given cluster is uniform in unsupervised scRNA-seq data analyses. Here, we proposed the concept of cluster purity and introduced a conceptually novel statistic, named ROGUE, to examine whether a given cluster is a pure cell population. (bioRxiv preprint link [here](https://www.biorxiv.org/content/10.1101/819581v1)).

## Installation Guide
**Installing dependency package**  
Before installing ROGUE, the “tidyverse” package should be installed first:
```
install.packages("tidyverse")
```
**Installing ROGUE**  
To install ROGUE, run:
```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("PaulingLiu/ROGUE")
```

## Tutorial
For more details and basic usage see following tutorials:
1.	[Guided Tutorial](https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html)

## Citation
If you use ROGUE in your research, please considering citing:
- [Liu et al., bioRxiv 2019](https://www.biorxiv.org/content/10.1101/819581v1)

## Contact
Please contact us:  
Baolin Liu: pauling.liu@pku.edu.cn  
Zemin Zhang: zemin@pku.edu.cn

## Copyright
©2019 Baolin Liu, Chenwei Li. [Zhang Lab](http://cancer-pku.cn/). All rights reserved.
