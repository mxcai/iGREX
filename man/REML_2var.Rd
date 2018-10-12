\name{REML_2var}
\alias{REML_2var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit variance components model (2 variance components) using MM algorithm
}
\description{Fit variance components model (2 variance components) using MM algorithm
}
\usage{
REML_2var(K, y, Z = NULL, maxIter = 1500, tol = 1e-06, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{Correlation (covariance) matrix.
}
  \item{y}{Response vector.
}
  \item{Z}{Covariates matrix. Default is a column of 1's.
}
  \item{maxIter}{Maximum iteration. Default is 1500.
}
  \item{tol}{Threshold to stop the EM iteration. Default is 1e-6
}
  \item{verbose}{Whether to show iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{List of outcomes:
  \item{beta0 }{Estimation of fixed effects.}
  \item{sb2 }{Estimation value of variance component.}
  \item{sy2 }{Estimation value of residual variance component.}
  \item{K }{Correlation (covariance) matrix.}
  \item{iter }{Total number of iteration used.}
  \item{covSig }{2 by 2 covariance matrix of 3 estimated variance components.}
  \item{h }{PVE estimate.}
  \item{se_h }{Standard error of PVE estimate.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Mingxuan Cai
}
\seealso{
  \code{\link{REML_3var}}, \code{\link{MoM_3var}}, \code{\link{MoM_2var}}
}
\examples{
}
