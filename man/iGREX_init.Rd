\name{iGREX_init}
\alias{iGREX_init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit iGREX stage one
}
\description{Fit iGREX stage one
}
\usage{
iGREX_init(prefix_eQTL_geno, prefix_GWAS, gene_expr, cov_eQTL, cov_GWAS, whCol, bw, subsample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prefix_eQTL_geno}{Prefix for eQTL genotype file with plink format (bim, bed).
}
  \item{prefix_GWAS}{Prefix for GWAS genotype and phenotype file with plink format (bim, bed, fam).
}
  \item{gene_expr}{Gene expression file with full name.
}
  \item{cov_eQTL}{Covariates file for eQTL data.
}
  \item{cov_GWAS}{Covariates file for GWAS data, e.g. top 10 PCs.
}
  \item{whCol}{Specify which phenotype is used in fam. For example, when whCol = 2, the seven-th column of fam file will be used as phenotype.
}
  \item{bw}{Number of downstream and upstream SNPs that are considered as cis-SNP of a gene.
}
  \item{subsample}{Number of subsample to use in order to ease the computation of Kinship matrix. Set to \0 if use the full sample.
}
}
\value{List of outcomes
  \item{idx_all }{Indices of cis-SNPs for each gene.}
  \item{idxinFile }{Indices of matched SNPs in GWAS file and eQTL file.}
  \item{ind }{Indicator vector takes value -1 or 1. -1 means the SNP is flipped.}
  \item{z }{Phenotype vector.}
  \item{covar1 }{Covariates of gene expression data (including a column of 1's).}
  \item{covar2 }{Covariates of GWAS data (including a vector of 1's).}
  \item{expr }{Matrix of all gene expressions. Each row is a gene, each column is a sample.}
  \item{snp_info }{Data frame containing the information of matched SNPs, e.g. Chr, rsname and BP.}
  \item{gene_info_all }{Data frame containing the information of all genes, e.g. lower, upper, genetype1, genetype2, TargetID, Chr.}
  \item{gene_info_match }{Data frame containing the information of matched genes, e.g. lower, upper, genetype1, genetype2, TargetID, Chr.}
  \item{out_param }{Matrix of PX-EM outputs for each gene. Columns are estimate of \sigma_{\beta}^2, estimate of \sigma_{y}^2, log-likelihood at convergence, number of cis_SNPs, gene level cis-heritability}
  \item{K }{GREX kinship matrix.}
  \item{K0 }{GREX kinship matrix without accounting for mearsurement error.}
  \item{subsample }{Subsample size.}
  \item{idx_subsample }{Indices of individuals selected if subsample is used.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Mingxuan Cai
}
\seealso{
\code{\link{iGREX}}, \code{\link{iGREX_Kg}}, \code{\link{lmm_pxem}}
}
\examples{
See Vegnette.
}
