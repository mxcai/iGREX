# library(MASS)
# Ghuber = function(u, k=30, deriv=0){
#   if(!deriv) return(pmin(1, k/abs(u)))
#   abs(u)<=k
# }

ldsc <- function(Z,r2, N, tau=NULL, W=NULL){
  if(is.null(W)) W = rep(1,length(Z))
  if(is.null(tau)) tau = (mean(Z^2)-1)/mean(N*r2)
  Wv = 1/(1+tau*N*r2)^2
  # id = which(Z^2>30)
  # if(length(id)>0) Wv[id] = sqrt(Wv[id])
  # rcf = as.vector( rlm(I(Z^2) ~ I(N*r2), weight=W*Wv, psi=Ghuber, k=30)$coef )
  # return(rcf)
  fit = lm(I(Z^2-1) ~ -1 + I(N*r2), weight=W*Wv)
  return(list(fit=summary(fit),Wv=Wv))
}

# ldsc1 <- function(Z,r2, N, tau=NULL, W=NULL){
#   if(is.null(W)) W = rep(1,length(Z))
#   if(is.null(tau)) tau = (mean(Z^2)-1)/mean(N*r2)
#   Wv = 1/(1+tau*N*r2)^2
#   # id = which(Z^2>30)
#   # if(length(id)>0) Wv[id] = sqrt(Wv[id])
#   # rcf = as.vector( rlm(I(Z^2) ~ I(N*r2), weight=W*Wv, psi=Ghuber, k=30)$coef )
#   # return(rcf)
#   fit = lm(I(Z^2) ~ I(N*r2), weight=W*Wv)
#   return(list(fit=summary(fit),W=Wv))
# }
