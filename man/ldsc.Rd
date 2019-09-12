\name{ldsc}
\alias{ldsc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ldsc for estimating heritability
}
\description{ldsc for estimating heritability using GWAS summary statistics
}
\usage{
ldsc(Z, r2, N, tau = NULL, W = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{A vector of z-scores.
}
  \item{r2}{Unbiased LD scores for analyzed SNPs. Should be of the same length as Z.
}
  \item{N}{Sample size for obtaining the summary statistics.
}
  \item{tau}{An initial guess of \eqn{h^2}. Default is \eqn{(mean(Z^2)-1)/mean(N*r2)}.
}
  \item{W}{Weighting for each SNP. Default is one.
}
}
\value{A list of outcomes
  \item{fit }{A summary of linear regression results. The heritability estimate is given by the intercept.}
  \item{Wv }{Regression weights, the reciprocal of conditional variance function \eqn{Var(\chi^2|l_j)}}
}
\references{
  Bulik-Sullivan, B. K., Loh, P. R., Finucane, H. K., Ripke, S., Yang, J., Patterson, N., ... & Schizophrenia Working Group of the Psychiatric Genomics Consortium. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nature genetics, 47(3), 291.
}

\author{
  Mingxuan Cai
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{MQS}}
}
