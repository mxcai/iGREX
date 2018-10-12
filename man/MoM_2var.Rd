\name{MoM_2var}
\alias{MoM_2var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit variance components model (2 variance components) using Method of Moments
}
\description{Fit variance components model (2 variance components) using Method of Moments
}
\usage{
MoM_2var(K, y, Z = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{Correlation (covariance) matrix.
}
  \item{y}{Response vector.
}
  \item{Z}{Covariates matrix. Default is a column of 1's.
}
}
\value{List of outcomes:
  \item{s12 }{Estimation value of variance component.}
  \item{s22 }{Estimation value of residual variance component.}
  \item{K }{Correlation (covariance) matrix.}
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
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{REML_2var}}, \code{\link{MoM_3var}}, \code{\link{MoM_3var}}
}
\examples{
}

