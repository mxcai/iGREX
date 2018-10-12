\name{lmm_pxem}
\alias{lmm_pxem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit linear mixed model by PX-EM
}
\description{
Fit linear mixed model by PX-EM
}
\usage{
lmm_pxem(y, X, W, tol, maxIter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Gene expression vector.
}
  \item{X}{Normalized genotype (cis-SNPs) matrix for eQTL.
}
  \item{W}{Covariates file for eQTL data.
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
  \code{\link{iGREX}}, \code{\link{iGREX_init}}, \code{\link{iGREX_Kg}}
}
\examples{
n <- 500
p <- 500
X <- matrix(rnorm(n*p),n,p)
beta <- rnorm(p,0,sqrt(0.5/p))

W <- matrix(1,n,1)
beta0 <- 0.1

y <- W\%*\%beta0 + X\%*\%beta + rnorm(n,0,sqrt(0.5))
fit <- lmm_pxem(y,X,W,1e-8,300)
}

