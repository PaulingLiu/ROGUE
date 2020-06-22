# ROGUE (Ratio of Global Unshifted Entropy)

## Contents

- [Overview](#overview)
- [Installation Guide](#installation-guide)
- [Tutorial](#tutorial)
- [Reproduction instructions](#Reproduction-instructions)
- [License](./LICENSE)
- [Citation](#citation)
- [Contact](#Contact)

## Overview
Often, it is not even clear whether a given cluster is uniform in unsupervised scRNA-seq data analyses. Here, we proposed the concept of cluster purity and introduced a conceptually novel statistic, named ROGUE, to examine whether a given cluster is a pure cell population.

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
1.	[Guided Tutorial](https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html) (It takes a few seconds to load the HTML file)

## Reproduction instructions
The scripts for producing all the quantitative results in our manuscript can be found in [scripts](./scripts).

## Citation
If you use ROGUE in your research, please considering citing:
- [Liu et al., Nature Communications 2020](https://www.nature.com/articles/s41467-020-16904-3)

## Contact
Please contact us:  
Baolin Liu: pauling.liu@pku.edu.cn  
Zemin Zhang: zemin@pku.edu.cn

## Copyright
©2019 Baolin Liu, Chenwei Li. [Zhang Lab](http://cancer-pku.cn/). All rights reserved.
