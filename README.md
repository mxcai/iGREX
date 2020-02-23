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
[The 'iGREX' vignette](https://github.com/mxcai/iGREX/blob/master/vignettes/iGREX.pdf?raw=true) provides a quick start for the usage of the package. The following help page also provides quick reference and examples:

```
library(iGREX)
package?iGREX
```

The `iGREXs` function requires only the GWAS summary statistics and the eQTL data for analysis, which greatly enhances the applicability of our software. The analysis in the vignette using GWAS summary statistics can be reproduced with open access data. The processed High-Density Lipoprotein cholesterol (HDL) GWAS summary statistics, Geuvadis gene expression data and 1000 Genomes genotype data used in the analysis can be downloaded from [this Dropbox link](https://www.dropbox.com/sh/xbq0a0or1nmcaef/AABcmzTxgWPJGpcCj9mYOzwma?dl=0). One can follow the instructions in the vignette to analyze these dataset.

Reproducibility
==========

All the simulation results can be reproduced by using the r scripts under IGREX-sim folder of thsi repository. To reproduce the results in our manuscript, one can first clone this repository and goes to the IGREX_sim folder:
```
git clone https://github.com/mxcai/iGREX
cd iGREX/IGREX_sim
mkdir results
```
Running the simulation.sh script will generate all the results for creating figures in the manuscrupt. Note that each chunk of the codes in the bash script takes hours to finish a single parameter setting. It is suggested to be run on the server.
```
chmod +x simulation.sh
./simulation.sh
```
One can also run only part of the codes at a time. For example, to produce the results shown in Figure 1e of the manuscript, use the following:
```
for rho in {0.1,0.3,0.5,0.8}
do
Rscript LD_rho.R 1000 4000 100 200 0.3 0.2 0.3 500 30 $rho
done
```
Once the results are obtained, one can use the r script generate_figures.R to produce figures in the manuscript.

References
==========

Mingxuan Cai, Lin S Chen, Jin Liu and Can Yang. IGREX for quantifying the impact of genetically regulated expression on phenotypes. NAR Genomics and Bioinformatics. https://doi.org/10.1093/nargab/lqaa010.


Development
==========

This R package is developed by Mingxuan Cai (mcaiad@ust.hk).
