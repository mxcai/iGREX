args <- commandArgs(trailingOnly = TRUE)

library(mvtnorm)
library(iGREX)
set.seed(10)

n1 <- as.numeric(args[1])
n2 <- as.numeric(args[2])
p1 <- as.numeric(args[3])   #number of SNPs in each gene
p2 <- as.numeric(args[4])   #number of genes
r <- 0     #number of overlapped snps in each gene
p <- (p1-r)*(p2-1)+p1
n_block <- 50
block_cut <- c(0,(1:n_block)*(p%/%n_block),p)     #block cutoff



sb2_true <- as.numeric(args[5])
sy2_true <- 1-sb2_true
sg2_true <- as.numeric(args[6])
sa2_true <- as.numeric(args[7])
f <- runif(p, 0.05, 0.5)
m <- as.numeric(args[8])

n_rep <- as.numeric(args[9])
out <- matrix(0,n_rep,4)

for(i in 1:n_rep){
  cat(i,"-th loop\n")
  X_tmp <- matrix(rnorm((n1+n2)*p),n1+n2,p)
  X <- matrix(1, n1+n2, p)
  X[t(t(X_tmp) - quantile(X_tmp, f^2)) < 0] <- 2
  X[t(t(X_tmp) - quantile(X_tmp, 1-(1-f)^2)) > 0] <- 0
  X <- scale(X)
  Y0 <- matrix(0,n1+n2,p2)

  for(g in 1:p2){
    beta <- rnorm(p1,0,sqrt(sb2_true/p1))
    Y0[,g] <- X[,((g-1)*(p1-r)+1):((g-1)*(p1-r)+p1)] %*% beta


  }
  Y <- Y0[1:n1,] + matrix(rnorm(n1*p2,0,sqrt(sy2_true)),n1,p2)

  X1 <- X[1:n1,]
  X2 <- X[(n1+1):(n1+n2),]

  t <- mean(diag(Y0[(n1+1):(n1+n2),]%*%t(Y0[(n1+1):(n1+n2),])))
  alpha <- as.matrix(rnorm(p2,0,sqrt(sg2_true/t)))
  z0 <- Y0[(n1+1):(n1+n2),] %*% alpha
  gamma <- as.matrix(rnorm(p,0,sqrt(sa2_true/p)))
  z1 <- X2%*%gamma
  z <- z0 + z1 + rnorm(n2,0,sqrt(var(z0+z1)))
  med_H_true <- var(z0)/var(z)

  z_score <- rep(0,p)
  for(j in 1:p){
    fit_lm <- lm(z~.,data = data.frame(z,X2[,j]))
    z_score[j] <- summary(fit_lm)$coefficients[2,3]
  }

  # Data: gene expr Y, phenotype z, genotype1 X1, genotype2 X2

  # fit LMM for step 1 and get K_g for each gene

  K <- K0 <- Km <- Km0 <- 0
  # K_diag <- vector("numeric",0)
  Km_diag <- vector("numeric",0)
  idx <- sample(1:n2,m,replace = F)
  q1_vec <- rep(0,p2)
  for(g in 1: p2){
    cat(g,"/",p2," gene\n")
    y_g <- Y[,g]
    X1tmp <- X1[,((g-1)*(p1-r)+1):((g-1)*(p1-r)+p1)]
    X2tmp <- X2[,((g-1)*(p1-r)+1):((g-1)*(p1-r)+p1)]
    ztmp <- z_score[((g-1)*(p1-r)+1):((g-1)*(p1-r)+p1)]
    W1 <- matrix(1,n1,1)
    W2 <- matrix(1,n2,1)
    fit_g <- iGREX_Kg(y_g,X1tmp,X2tmp,W1,1e-5,500)
    K <- K + fit_g$K_g
    K0 <- K0 + fit_g$K_g0

    q1_vec[g] <- t(ztmp/sqrt(n2))%*%fit_g$weight%*%ztmp/sqrt(n2) / p1

    fitrd_g <- iGREX_Kg(y_g,X1tmp,X2tmp[idx,],W1,1e-5,500)
    Km <- Km + fitrd_g$K_g
    Km0 <- Km0 + fitrd_g$K_g0
    Km_diag <- c(Km_diag,sum(diag(fitrd_g$K_g)))

  }
  q2_vec <- (z_score/sqrt(n2))^2

  mdiag <- mean(diag(K))
  K <- K/mdiag

  mdiagm <- mean(diag(Km))
  Km <- Km/mdiagm

  # X2s <- scale(X2)
  Ka <- X2 %*% t(X2) / ncol(X2)
  Xm <- scale(X2[idx,])
  Kma <- Xm %*% t(Xm) / ncol(Xm)

  # MoM using summary statisitcs
  trK1 <- sum(diag(Km))
  trK2 <- sum(diag(Kma))

  trK12 <- sum(Km^2)
  trK22 <- sum(Kma^2)
  trK1K2 <- sum(Km*Kma)

  c <- 1
  S <- matrix(0,2,2)
  S[1,1] <- (trK12-trK1^2/(m-c))/(m-c)^2
  S[1,2] <- S[2,1] <- (trK1K2-trK1*trK2/(m-c))/(m-c)^2
  S[2,2] <- (trK22-trK2^2/(m-c))/(m-c)^2
  q_ss <- c(sum(q1_vec)/mdiagm - 1/n2, sum(q2_vec)/p - 1/n2)
  invS <- solve(S)
  med_H_ss <- invS%*%q_ss

  group1 <- rep(0,p2)
  group2 <- rep(0,p)
  idx_group <- 1
  n_block <- length(block_cut)-1

  for(j in 1:n_block){
    tmp1 <- rep(FALSE,p2)
    tmp2 <- rep(FALSE,p)
    for(g in 1: p2){
      # cat(g,"/",p2," gene\n")
      gstart <- (g-1)*(p1-r)+1
      gend <- (g-1)*(p1-r)+p1
      if(gstart>=(block_cut[j]+1) & gend<=block_cut[j+1]){
        tmp1[g] <- TRUE
      }
    }
    tmp2[(block_cut[j]+1):block_cut[j+1]] <- TRUE

    if(sum(tmp1!=0)&sum(tmp2!=0)){
      group1[tmp1] <- idx_group
      group2[tmp2] <- idx_group
      idx_group <- idx_group+1
    }
  }
  ngroup <- idx_group-1
  qj <- sapply(1:ngroup,function(j){
    tmp1 <- group1==j
    tmp2 <- group2==j

    q1 <- sum(q1_vec[tmp1])
    q2 <- sum(q2_vec[tmp2])
    c(q1,q2,sum(Km_diag[tmp1])/m,sum(tmp2))
  })
  t1 <- sum(Km_diag[group1!=0])/m
  pp <- sum(group2!=0)

  q_j <- (c(sum(q1_vec[group1!=0]),sum(q2_vec[group2!=0])) - qj[1:2,])/(c(t1,pp)-qj[3:4,]) - 1/n2

  var_h <- invS %*% var(t(q_j)) %*% invS * (ngroup-1)


  out[i,1:2] <- med_H_ss
  out[i,3:4] <- sqrt(diag(var_h))
}

out <- data.frame(out)
names(out) <- c("PVEg_ss", "PVEa_ss","se_PVEg_ss", "se_PVEa_ss")
write.table(out,file=paste("./results/MoM_subsample",m,"_prop",sg2_true*10,sa2_true*10,"_SNRy",sb2_true,"_n",n1,"_",n2,".txt",sep=""),quote = F,col.names = T,row.names = F)
