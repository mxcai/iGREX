\name{MoM_3var_ss}
\alias{MoM_3var_ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit iGREXs stage two
}
\description{Fit iGREXs stage two
}
\usage{
MoM_3var_ss(obj, K2, sd_method = "LD_block")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{The object of iGREX stage one result. The outout of \code{\link{iGREXs_init}}.
}
  \item{K2}{Kinship matrix (GRM).
}
  \item{sd_method}{Method used to estimate the standard error. Currently support LD bloks only (default).
}
}
\value{List of outcomes: return(list(H=H,S=S,q=q,cov_h=var_h))
  \item{PVE }{2 by 3 matrix of PVE estimates and their standard errors. Columns correspond to \eqn{PVE_{GREX}}, \eqn{PVE_{Alternative}} and \eqn{PVE_{GREX}/h^2_t}.}
  \item{S }{Coefficiants of the normal equation.}
  \item{q }{RHS of the normal equation.}
  \item{cov_PVE }{2 by 2 covariance matrix of the \eqn{PVE_{GREX}} and \eqn{PVE_{Alternative}} estimates.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Mingxuan Cai
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{iGREXs}}, \code{\link{iGREXs_init}}, \code{\link{MoM_3var}}
}
