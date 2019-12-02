args <- commandArgs(trailingOnly = TRUE)

library(mvtnorm)
library(iGREX)
set.seed(10)

n1 <- as.numeric(args[1])
n2 <- as.numeric(args[2])
p1 <- as.numeric(args[3])   #number of SNPs in each gene
p2 <- as.numeric(args[4])   #number of genes
p <- p1*p2

sb2_true <- as.numeric(args[5])
sy2_true <- 1-sb2_true
sg2_true <- as.numeric(args[6])
sz2_true <- as.numeric(args[7])

m <- as.numeric(args[8])

n_rep <- as.numeric(args[9])
out <- matrix(0,n_rep,5)

for(i in 1:n_rep){
  cat(i,"-th loop\n")

  X <- matrix(rnorm((n1+n2)*p1*p2),n1+n2,p1*p2)
  X <- scale(X)
  Y0 <- matrix(0,n1+n2,p2)

  for(g in 1:p2){
    beta <- rnorm(p1,0,sqrt(sb2_true/p1))
    Y0[,g] <- X[,(g*p1-p1+1):(g*p1)] %*% beta


  }
  Y <- Y0[1:n1,] + matrix(rnorm(n1*p2,0,sqrt(sy2_true)),n1,p2)

  X1 <- X[1:n1,]
  X2 <- X[(n1+1):(n1+n2),]

  t <- mean(diag(Y0[(n1+1):(n1+n2),]%*%t(Y0[(n1+1):(n1+n2),])))
  alpha <- as.matrix(rnorm(p2,0,sqrt(sg2_true/t)))
  z0 <- Y0[(n1+1):(n1+n2),] %*% alpha


  z <- z0 + rnorm(n2,0,sqrt(sz2_true))
  med_H_true <- var(z0)/var(z)

  z_score <- rep(0,p)
  for(j in 1:p){
    fit_lm <- lm(z~.,data = data.frame(z,X2[,j]))
    z_score[j] <- summary(fit_lm)$coefficients[2,3]
  }

  # Data: gene expr Y, phenotype z, genotype1 X1, genotype2 X2

  # fit LMM for step 1 and get K_g gor each gene

  K <- K0 <- Km <- Km0 <- 0
  idx <- sample(1:n2,m,replace = F)
  q1_vec <- rep(0,p2)
  q0_vec <- rep(0,p2)
  z_TWAS <- rep(0,p2)
  Y_hat <- matrix(0,n2,p2)
  for(g in 1: p2){
    cat(g,"/",p2," gene\n")
    y_g <- Y[,g]
    X1tmp <- X1[,(g*p1-p1+1):(g*p1)]
    X2tmp <- X2[,(g*p1-p1+1):(g*p1)]
    ztmp <- z_score[(g*p1-p1+1):(g*p1)]
    W1 <- matrix(1,n1,1)
    W2 <- matrix(1,n2,1)
    fit_g <- iGREX_Kg(y_g,X1tmp,X2tmp,W1,1e-5,500)
    K <- K + fit_g$K_g
    K0 <- K0 + fit_g$K_g0

    q1_vec[g] <- t(ztmp/sqrt(n2))%*%fit_g$weight%*%ztmp/sqrt(n2) / p1
    q0_vec[g] <- sum(ztmp/sqrt(n2)*fit_g$mub)^2 / p1

    z_TWAS[g] <- sum(ztmp*fit_g$mub)/sqrt(t(fit_g$mub)%*%cor(X2tmp[idx,])%*%fit_g$mub)
    Y_hat[,g] <- X2tmp[idx,]%*%fit_g$mub

    fitrd_g <- iGREX_Kg(y_g,X1tmp,X2tmp[idx,],W1,1e-5,500)
    Km <- Km + fitrd_g$K_g
    Km0 <- Km0 + fitrd_g$K_g0

  }


  mdiag <- mean(diag(K))
  K <- K/mdiag

  mdiagm <- mean(diag(Km))
  Km <- Km/mdiagm

  mdiagm0 <- mean(diag(Km0))
  Km0 <- Km0/mdiagm0

  # REML
  REML <- REML_2var(K,z)

  # exact estimate by MoM
  trK <- sum(diag(K))
  tr2K <- trK^2
  trK2 <- sum(K^2)

  temp <- K - trK/(n2-1)*diag(n2)
  denom <- trK2 - tr2K/(n2-1)
  sg2 <- drop((t(z) %*% temp %*% z) / denom)
  sz2 <- drop((t(z) %*% (trK2/(n2-1)*diag(n2) - K*trK/(n2-1)) %*% z) / denom)

  MoM_H <- trK*sg2/sum(z^2)

  # MoM using summary statisitcs
  trK_ss <- sum(diag(Km))
  tr2K_ss <- trK_ss^2
  trK2_ss <- sum(Km^2)

  S <- (trK2_ss- tr2K_ss/(m-1))/(m-1)^2
  q_ss <- sum(q1_vec)/mdiagm - 1/n2
  MoM_H_ss <- q_ss/S

  # ldsc-RHOGEn
  corY2 <- cor(Y_hat)^2
  r2_unbiased <- colSums(corY2-(1-corY2)/m)
  fit_ldsc <- ldsc(z_TWAS,r2_unbiased/p2,n2)

  # MoM using summary statisitcs without accounting for uncertainty

  trK0_ss <- sum(diag(Km0))
  tr2K0_ss <- trK0_ss^2
  trK02_ss <- sum(Km0^2)

  S <- (trK02_ss- tr2K0_ss/(m-1))/(m-1)^2
  q_ss <- sum(q0_vec)/mdiagm0 - 1/n2
  MoM_H_ss0 <- q_ss/S

  out[i,] <- c(REML=REML$h,MoM_H,MoM_H_ss,MoM_H_ss0,fit_ldsc$fit$coefficients[1])
}

colnames(out) <- c("REML","MoM","iGREX_ss","iGREX_ss0","RHOGE")
write.table(out,file=paste("./results/iGREXvsRHOGE_SNRz",sg2_true,"_SNRy",sb2_true,"_n",n1,"_",n2,".txt",sep=""),quote = F,col.names = T,row.names = F)

