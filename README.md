iGREX
===

Quantifying the impact of genetically regulated expression on complex traits and diseases.

Installation
===========

To install the development version of bivas, it's easiest to use the 'devtools' package. Note that bivas depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("mxcai/iGREX")
```

Usage
===========
[The 'iGREX' vignette](https://github.com/mxcai/iGREX/blob/master/vignettes/iGREX_real_data.pdf) provides a quick start for the usage of the package. The following help page also provides quick reference and examples:

```
library(iGREX)
package?iGREX
```

Reproducibility
==========

All the simulation results can be reproduced by using the code at [sim-iGREX](https://github.com/mxcai/sim-iGREX).
The `iGREXs` function requires only the GWAS summary statistics and the eQTL data for analysis, which greatly enhances the applicability of our software. The analysis in the vignette using GWAS summary statistics can be reproduced using open access data. The processed High-Density Lipoprotein cholesterol (HDL) GWAS summary statistics, Geuvadis gene expression data and 1000 Genomes genotype data used in the analysis can be downloaded from [this Dropbox link](https://www.dropbox.com/sh/xbq0a0or1nmcaef/AABcmzTxgWPJGpcCj9mYOzwma?dl=0).

References
==========

Mingxuan Cai, Lin S Chen, Jin Liu and Can Yang. IGREX for quantifying the impact of genetically regulated expression on phenotypes. biorxiv: https://www.biorxiv.org/content/10.1101/546580v2.


Development
==========

This R package is developed by Mingxuan Cai (mcaiad@ust.hk).
