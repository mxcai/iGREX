\name{iGREXs}
\alias{iGREXs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit iGREX-s using GWAS summary statistics
}
\description{Fit iGREX-s using GWAS summary statistics
}
\usage{
iGREXs(prefix_eQTL_geno, prefix_GWAS, gene_expr, Z_score, cov_eQTL = "", cov_GWAS = "", trans_eQTL = "", Ka, whCol = 1, bw = 5e+05, sd_method = "LD_block")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prefix_eQTL_geno}{Prefix for eQTL genotype file with plink format (bim, bed, fam).
}
  \item{prefix_GWAS}{Prefix for GWAS genotype and phenotype file with plink format (bim, bed, fam).
}
  \item{gene_expr}{Gene expression file with full name.
}
  \item{Z_score}{Z-score file of GWAS study.
}
  \item{cov_eQTL}{Covariates file for eQTL data.
}
  \item{cov_GWAS}{Covariates file for GWAS data, e.g. top 10 PCs.
}
  \item{trans_eQTL}{File containg the trans-eQTLs to be included in the study.
}
  \item{Ka}{Kinship matrix (GRM).
}s
  \item{bw}{Number of downstream and upstream SNPs that are considered as cis-SNP of a gene.
}
  \item{sd_method}{Method used to estimate the standard error. Currently support LD bloks only (default).
}
}

\value{List of outcomes:
  \item{fit_init }{List of Stage one results. Same as the output of \code{\link{iGREXs_init}}.}
  \item{output }{List of Stage two results. Same as the output of \code{\link{MoM_3var_ss}}.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Mingxuan Cai
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{MoM_3var_ss}}, \code{\link{iGREXs_init}}
}
\examples{
See Vegnette.
}
