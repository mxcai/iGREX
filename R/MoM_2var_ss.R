MoM_2var_ss <- function(Km,zs,n,Z=NULL,group=NULL){
  m <- nrow(Km)
  p <- length(zs)
  ngroup <- length(unique(group))

  if(is.null(Z)){
    Z <- matrix(1,m,1)
    M <- diag(m) - matrix(1/m,m,m)
    MK <- Km
  } else{
    Z <- cbind(1,Z)
    M <- diag(m) - Z%*%solve(t(Z)%*%Z)%*%t(Z)
    MK <- M%*%Km
  }
  q <- ncol(Z)

  trKm <- sum(diag(MK))
  tr2Km <- trKm^2
  trKm2 <- sum(MK^2)


  S <- (trKm2-tr2Km/(m-q)) / (m-q)^2
  # zz <- sum((zs/sqrt(n))^2)
  # c <- zz/p - 1/(median(n))
  zz <- zs^2/n - 1/n
  c <- sum(zz) / p

  h <- c/S

  if(is.null(group)){
    # compute standard error by Jackknife
    # zj <- (zs/sqrt(n))^2
    # c_jf <- (zz-zj)/(p-1) - 1/(median(n))
    c_jf <- (sum(zz)-zz)/(p-1)


    var_h <- var(c_jf)/S/S*(p-1)
  } else{
    # zj <- sapply(1:ngroup,function(j){
    #
    #   tmp <- sum((zs[group==j]/sqrt(n[group==j]))^2)
    #   c(tmp,sum(group==j))
    # })
    # c_jf <- (zz-zj[1,])/(p-zj[2,]) - 1/(median(n))
    #
    # var_h <- var(c_jf)/S/S*(ngroup-1)
    zj <- sapply(1:ngroup,function(j){

      tmp <- sum(zz[group==j])
      c(tmp,sum(group==j))
    })
    c_jf <- (sum(zz)-zj[1,])/(p-zj[2,])

    var_h <- var(c_jf)/S/S*(ngroup-1)
  }


  # ret <- list(h=h,Km=Km,se_h=sqrt(var_h),c_jf=c_jf,S=S,c=c)
  ret <- list(h=h,Km=Km,se_h=sqrt(var_h))
}
