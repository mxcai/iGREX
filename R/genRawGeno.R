library(mvtnorm)
genRawGeno <- function(maf, p, G, rho, n){
  sigma <- rho^(abs(outer(1:p, 1:p, "-")))

  X_tmp <- replicate(G,rmvnorm(n,mean=rep(0,p),sigma=sigma,method = "chol"))
  dim(X_tmp) <- c(n,p*G)
  X <- matrix(1, n, p*G)
  X[t(t(X_tmp) - mapply(quantile,x=data.frame(X_tmp),probs=maf^2)) < 0] <- 2
  X[t(t(X_tmp) - mapply(quantile,x=data.frame(X_tmp),probs=1-(1-maf)^2)) > 0] <- 0
  return(X)
}
