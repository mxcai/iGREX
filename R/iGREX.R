iGREX <- function(prefix_eQTL_geno, prefix_GWAS, gene_expr, cov_eQTL="", cov_GWAS="", trans_eQTL="", Ka, whCol=1, bw=500000, subsample=0, method="REML"){
  cat("########## Stage One starts ##########\n")
  fit_init <- iGREX_init(prefix_eQTL_geno,prefix_GWAS,gene_expr,cov_eQTL,cov_GWAS,trans_eQTL,whCol,bw,subsample)

  cat("########## Stage Two starts ##########\n")
  G <- nrow(fit_init$gene_info_match)
  if(method=="REML"){
    fit_iGREX <- REML_3var(fit_init$K/G,Ka,fit_init$z,fit_init$cov_GWAS[,-1])
  } else if(method=="MoM"){
    fit_iGREX <- MoM_3var(fit_init$K/G,Ka,fit_init$z,fit_init$cov_GWAS[,-1])
  }
  return(list(fit_init=fit_init,output=fit_iGREX))
}


iGREXs <- function(prefix_eQTL_geno, prefix_GWAS, gene_expr, Z_score, cov_eQTL="", cov_GWAS="", trans_eQTL="", Ka, bw=500000, sd_method="LD_block"){
  cat("########## Stage One starts ##########\n")
  fit_init <- iGREXs_init(prefix_eQTL_geno, prefix_GWAS, gene_expr, Z_score, cov_eQTL, cov_GWAS, trans_eQTL, bw)

  cat("########## Stage Two starts ##########\n")
  fit_iGREX <- MoM_3var_ss(fit_init,Ka,sd_method = sd_method)

  return(list(fit_init=fit_init,output=fit_iGREX))
}
