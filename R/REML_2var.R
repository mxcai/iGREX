# ridge regression implemented by MM
REML_2var <- function(K,y,Z=NULL,maxIter=1500,tol=1e-6,verbose=T){
  # p <- ncol(X)
  n <- length(y)

  if(is.null(Z)){
    Z <- matrix(1,n,1)
  } else{
    Z <- cbind(1,Z)
  }

  #means of X, y and z
  # Xm <- colMeans(X)
  ym <- mean(y)


  #calculate sth in advance
  if(ncol(Z)==1){
    SZy <- ym
    # SZX <- Xm
    # X <- scale(X,scale = F,center = T)
    # y <- y - ym
  } else {
    ZZ <- t(Z)%*%Z
    SZy <- solve(ZZ,t(Z)%*%y)
    # SZX <- solve(ZZ,t(Z)%*%X)
    # X <- X - Z%*%SZX
    # y <- y - Z%*%SZy
  }


  # XX <- X%*%t(X)
  # Xy <- t(X)%*%y

  eigenK <- eigen(K)
  eVal <- eigenK$values
  eVec <- eigenK$vectors

  U <- eVec

  yt <- t(U)%*%y
  Zt <- t(U)%*%Z

  #initialize
  sy2 <- drop(var(y))
  sb2 <- drop(var(y)) #/ p
  # mu <- matrix(0,p,1)
  # beta0 <- SZy - SZX %*% mu
  lb0 <- rep(0,maxIter)
  lb0[1] <- -Inf

  for(iter in 2:maxIter){

    d <- 1/(sb2*eVal + sy2)
    beta0 <- solve(t(Zt)%*%diag(d)%*%Zt) %*% t(Zt) %*% diag(d) %*% yt
    res <- yt - Zt %*% beta0

    sb2 <- sb2 * sqrt(sum(res^2 * d^2 * eVal) / sum(d * eVal))
    sy2 <- sy2 * sqrt(sum(res^2 * d^2) / sum(d))

    #evaluate lower bound
    lb0[iter] <- get_lb_MM(d,res,beta0)

    if(verbose){
      cat(iter,"-th iteration, lower bound = ",lb0[iter]," ,diff=",lb0[iter]-lb0[iter-1],",sb2=",sb2,",sy2=",sy2,"\n")
    }
    if(abs(lb0[iter]-lb0[iter-1])<tol){
      lb0 <- lb0[1:iter]
      break
    }
  }

  invSigy <- solve(sb2*K+sy2*diag(n))
  invSigyK <- invSigy%*%K
  FIM <- matrix(0,2,2)
  FIM[1,1] <- sum(invSigyK^2) / 2
  FIM[2,2] <- sum(invSigy^2) / 2
  FIM[1,2] <- FIM[2,1] <- sum(invSigyK*invSigy) / 2
  covSig <- solve(FIM)  #inverse of FIM

  # calculate h = sb2*trK / (sb2*trK + sy2*n)
  trK <- sum(diag(K))
  var_total <- sb2*trK + sy2*n
  h <- sb2*trK / var_total
  gh <- c(sy2*trK*n/var_total^2, -sb2*trK*n/var_total^2)
  se_h <- sqrt(t(gh) %*% covSig %*% gh)

  bayesRegMM <- list(beta0=beta0,sb2=sb2,sy2=sy2,K=K,iter=iter,covSig=covSig,h=h,se_h=se_h)


  attr(bayesRegMM,"class") <- "bayesRegMM"
  bayesRegMM
}


####################################################################################################
get_lb_MM <- function(d,res,beta0){
  llh <- sum(log(d))/2 - sum(res^2*d)/2
}
