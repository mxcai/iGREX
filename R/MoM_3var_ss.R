MoM_3var_ss <- function(obj,K2,sd_method="LD_block"){
  m <- nrow(K2)
  n_gene <- length(obj$q1_vec)
  n_SNP <- length(obj$q2_vec)

  K1 <- obj$K
  K_diag <- obj$K_diag
  q <- obj$q
  q1_vec <- obj$q1_vec
  q2_vec <- obj$q2_vec
  Z <- obj$covar2
  idx_all <- obj$idx_all

  if(ncol(Z)==1){
    # Z <- matrix(1,m,1)
    M <- diag(m) - matrix(1/m,m,m)
    MK1 <- K1
    MK2 <- K2
  } else{
    # Z <- cbind(1,Z)
    M <- diag(m) - Z%*%solve(t(Z)%*%Z)%*%t(Z)
    MK1 <- M%*%K1
    MK2 <- M%*%K2
  }
  c <- ncol(Z)

  trK1 <- sum(diag(MK1))
  trK2 <- sum(diag(MK2))

  trK12 <- sum(MK1^2)
  trK22 <- sum(MK2^2)
  trK1K2 <- sum(MK1*MK2)

  S <- matrix(0,2,2)
  S[1,1] <- (trK12-trK1^2/(m-c))/(m-c)^2
  S[1,2] <- S[2,1] <- (trK1K2-trK1*trK2/(m-c))/(m-c)^2
  S[2,2] <- (trK22-trK2^2/(m-c))/(m-c)^2

  invS <- solve(S)
  h <- invS%*%q

  if(sd_method=="LD_block"){
    cat("Assigning SNPs to LD Blocks...\n")
    block <- read.table(system.file("extdata", "fourier_ls-all.bed", package = "medH"),header = T)
    group_gene <- rep(0,n_gene)
    group_SNP <- rep(0,n_SNP)
    snp_info <- obj$snp_info[unique(obj$idx_all[,1]),]
    idx_group <- 1
    for(i in 1:22){
      block_i <- block[block$chr==paste("chr",i,sep=""),]
      # chr_i <- obj$gene_info_match[obj$gene_info_match$Chr==i,]
      n_block <- nrow(block_i)    # number of blocks within each CHR

      for(j in 1:n_block){
        # genes entirely within j-th block in i-th CHR
        tmp1 <- with(obj$gene_info_match,Chr==i & (lower-500000)>=block_i$start[j] & (upper+500000)<block_i$stop[j])
        # SNPs within j-th block in i-th CHR
        tmp2 <- with(snp_info,Chr==i & BP>=block_i$start[j] & BP<block_i$stop[j])

        # the block is used in calculating sd if it's non-empty
        if(sum(tmp1!=0)&sum(tmp2!=0)){
          group_gene[tmp1] <- idx_group
          group_SNP[tmp2] <- idx_group
          idx_group <- idx_group+1
        }
        # if(sum(tmp2!=0)){
        #   group_SNP[tmp2] <- idx_group
        #   idx_group <- idx_group+1
        # }
      }
    }
    ngroup <- idx_group-1
    qj <- sapply(1:ngroup,function(j){
      tmp1 <- group_gene==j
      tmp2 <- group_SNP==j

      q1 <- sum(q1_vec[tmp1])
      q2 <- sum(q2_vec[tmp2])
      c(q1,q2,sum(K_diag[tmp1])/m,sum(tmp2))
    })
    t1 <- sum(K_diag[group_gene!=0])/m
    p1 <- sum(group_SNP!=0)

    # qj[3,] <- qj[3,]/n_gene1*obj$mdiag
    q_j <- (c(sum(q1_vec[group_gene!=0]),sum(q2_vec[group_SNP!=0])) - qj[1:2,])/(c(t1,p1)-qj[3:4,]) - 1/median(obj$z_info_match$N)

    var_h <- invS %*% var(t(q_j)) %*% invS * (ngroup-1)

    # calculate se of h_med = h1 / (h1 + h2) using Delta method
    gh <- c(h[2]/(h[1] + h[2])^2, -h[1]/(h[1] + h[2])^2)
    se_h_med <- sqrt(t(gh) %*% var_h %*% gh)
  }

  H <- matrix(c(h,h[1]/sum(h),sqrt(diag(var_h)),se_h_med),2,3,byrow = T)
  colnames(H) <- c("PVEg","PVEa","prop")
  rownames(H) <- c("heritability","se")

  return(list(PVE=H,S=S,q=q,cov_PVE=var_h))
}
