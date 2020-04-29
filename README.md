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

Application of IGREX-s to GTEx data
===========
To apply IGREX-s to GTEx data with GWAS summary statistics, you need the following input files:
1. The genotype file of GTEx samples, which should be applied through dbGap. Find instructtions in https://gtexportal.org/home/protectedDataAccess. In the function iGREXs, argument `prefix_eQTL_geno` takes the prefix of the genotype file in plink format.
2. The expression file of interested tissue, preprocessed and normalized. You can download this data from the GTEx websit. The v7 release used in our paper can be found at https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_expression_matrices.tar.gz. In the function iGREXs, argument `gene_expr` takes the full name of the expression file. File format should follow the one in the vignette.
3. The covarites file of eQTL, including genotype principal components and PEER factors as provided in the GTEx website. The v7 release used in our paper can be found at https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz. In the function iGREXs, argument `cov_eQTL` takes the full name of the eQTL covariates file.
4. GWAS summary statistics in the required format as shown in the vignette. In the function iGREXs, argument `Z_score` takes the full name of the eQTL covariates file.
5. A genotype file as LD reference for the GWAS. In practice you can simply use the same genotype file of GTEx samples in 1. In the function iGREXs, argument `prefix_GWAS` takes the prefix of the genotype file in plink format.
6. The covarites file of the GWAS LD reference, age, sex, genotype PC's, etc. If you are using the genotypes of GTEx samples as LD reference, This file contains the PC columns of eQTL covarites in 3. 
7. A kinship matrix of the LD reference.

Suppose we are using GTEx genotypes as GWAS LD reference, IGREX-s can be fitted by:

```
file_geno_gtex = "GTEx_qc_hm3" # prefix of GTEx genotype plink files (GTEx_qc_hm3.bim/GTEx_qc_hm3.bed/GTEx_qc_hm3.fam)
file_geno_GWASref = "GTEx_qc_hm3" # use GTEx genotypes as GWAS LD reference
file_expr = "Liver_gene_expression.txt" # full name of normalized gene expression matrix file 
file_covar_gtex = "Liver_cov.txt" # covariates of GTEx eQTL data
file_cov_GWASref = "GTEx_qc_hm3_pc5.eigenvec" # covariates of GWAS LD reference 
file_z = "HDL_summary.txt" # summary statistics file

fit_IGREXs <- iGREXs(prefix_eQTL_geno=file_qtlgeno,prefix_GWAS=file_GWASgeno,gene_expr=file_expression,Z_score=file_z,
cov_eQTL=file_qtlcov,cov_GWAS=file_GWAScov,Ka=K,bw = 500000)
```


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
