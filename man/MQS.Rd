\name{MQS}
\alias{MQS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MQS for estimating heritability
}
\description{MQS for estimating heritability using GWAS summary statistics
}
\usage{
MQS(file_z, file_ref, file_cov = NULL, K = NULL, X = NULL, sd_method = "LD_block")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_z}{Summary statistics file.
}
  \item{file_ref}{Prefix of reference Plink files.
}
  \item{file_cov}{Covariates file.
}
  \item{K}{Kinship matrix. (Optional. If not provided, the kinship matrix will be computed from the reference samples in the Plink files.)
}
  \item{X}{Genotype matrix. (Optional. If not provided, X will be read from the reference samples in the Plink files, and the kinship matrix is computed from X.)
}
  \item{sd_method}{String specifies the method for computing standard error. One of "LD_block", "Chromosome" or "Jackknife". Default is "LD_block".
}
}

\value{List of outcomes
  \item{h }{heritability estimation.}
  \item{K }{Kinship matrix.}
}
\references{
  Zhou, X. (2017). A unified framework for variance component estimation with summary statistics in genome-wide association studies. The annals of applied statistics, 11(4), 2027.
}

\author{
  Mingxuan Cai
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ldsc}}
}

