\name{iGREX}
\alias{iGREX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit iGREX
}
\description{Fit iGREX
}
\usage{
iGREX(prefix_eQTL_geno, prefix_GWAS, gene_expr, cov_eQTL, cov_GWAS, Ka, whCol = 1, bw = 5e+05, subsample = 0, method = "REML")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prefix_eQTL_geno}{Prefix for eQTL genotype file with plink format (bim, bed, fam).
}
  \item{prefix_GWAS}{Prefix for GWAS genotype and phenotype file with plink format (bim, bed, fam).
}
  \item{gene_expr}{Gene expression file with full name.
}
  \item{cov_eQTL}{Covariates file for eQTL data.
}
  \item{cov_GWAS}{Covariates file for GWAS data, e.g. top 10 PCs.
}
  \item{Ka}{Kinship matrix (GRM).
}
  \item{whCol}{Specify which phenotype is used in fam. For example, when whCol = 2, the seven-th column of fam file will be used as phenotype.
}
  \item{bw}{Number of downstream and upstream SNPs that are considered as cis-SNP of a gene.
}
  \item{subsample}{Number of subsample to use in order to ease the computation of Kinship matrix. Set to \0 if use the full sample.
}
  \item{method}{String specifies the method used in stage two, should be either "MoM" or "REML. Default is "REML".
}
}
\value{List of outcomes:
  \item{fit_init }{List of Stage one results. Same as the output of \code{\link{iGREX_init}}.}
  \item{output }{List of Stage two results. Same as the output of \code{\link{REML_3var}} or \code{\link{REML_3var}}.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Mingxuan Cai
}
\seealso{
 code{\link{REML_3var}}, code{\link{MoM_3var}}, code{\link{iGREX_init}}
}
\examples{
See Vegnette.
}
