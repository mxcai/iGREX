\name{iGREX_Kg}
\alias{iGREX_Kg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate GREX kinship matrix for one gene
}
\description{Calculate GREX kinship matrix
}
\usage{
iGREX_Kg(y, X1, X2, W1, tol, maxIter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Gene expression vector.
}
  \item{X1}{Normalized genotype (cis-SNPs) matrix for eQTL.
}
  \item{X2}{Normalized genotype (cis-SNPs) matrix for GWAS.
}
  \item{W1}{Covariates file for eQTL data.
}
  \item{tol}{Threshold to stop the EM iteration.
}
  \item{maxIter}{Maximum iteration.
}
}
\value{List of outcomes:
  \item{beta0 }{Estimation of fixed effects.}
  \item{sigma2y }{Estimation of residual variance component.}
  \item{sigma2beta }{Estimation of genetic variance component.}
  \item{Sigb }{Posterior covariance matrix of beta.}
  \item{K_g }{GREX kinship matrix of one gene.}
  \item{K_g0 }{GREX kinship matrix of one gene without accounting for mearsurement error.}
  \item{mub }{Posterior mean vector of beta.}
  \item{loglik }{Log-likelihood value at convrgence.}
  \item{iter }{Total number of iteration used.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Mingxuan Cai
}
\seealso{
\code{\link{iGREX}}, \code{\link{iGREX_init}}, \code{\link{lmm_pxem}}
}
\examples{
See vegnette.
}
