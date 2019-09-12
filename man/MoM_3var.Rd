\name{MoM_3var}
\alias{MoM_3var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit variance components model (3 variance components) using Method of Moments
}
\description{Fit variance components model (3 variance components) using Method of Moments
}
\usage{
MoM_3var(K1, K2, y, Z = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K1}{Correlation (covariance) matrix 1.
}
  \item{K2}{Correlation (covariance) matrix 2.
}
  \item{y}{Response vector.
}
  \item{Z}{Covariates matrix. Default is a column of 1's.
}
}
\value{List of outcomes:
  \item{s12 }{Estimation value of variance component 1.}
  \item{s22 }{Estimation value of variance component 2.}
  \item{se2 }{Estimation value of residual variance component.}
  \item{K1 }{Correlation (covariance) matrix 1.}
  \item{K2 }{Correlation (covariance) matrix 2.}
  \item{covSig }{3 by 3 covariance matrix of 3 estimated variance components.}
  \item{PVE }{2 by 4 matrix of PVE estimates and their standard errors. Columns correspond to \eqn{PVE_{GREX}}, \eqn{PVE_{Alternative}}, \eqn{h^2_t} and \eqn{PVE_{GREX}/h^2_t}.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Mingxuan Cai
}
\seealso{
\code{\link{iGREX}}, \code{\link{iGREX_Kg}}, \code{\link{iGREX_init}}, \code{\link{REML_3var}}, \code{\link{REML_2var}}, \code{\link{MoM_2var}}
}
\examples{
See vegnette.
}
