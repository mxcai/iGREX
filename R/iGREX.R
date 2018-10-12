iGREX <- function(prefix_eQTL_geno, prefix_GWAS, gene_expr, cov_eQTL,cov_GWAS, Ka, whCol=1, bw=500000, subsample=0, method="REML"){
  cat("########## Stage One starts ##########\n")
  fit_init <- iGREX_init(prefix_eQTL_geno,prefix_GWAS,gene_expr,cov_eQTL,cov_GWAS,whCol,bw,subsample)

  cat("########## Stage Two starts ##########\n")
  G <- nrow(fit_init$gene_info_match)
  if(method=="REML"){
    fit_iGREX <- REML_3var(fit_init$K/G,Ka,fit_init$z,fit_init$cov_GWAS[,-1])
  } else if(method=="MoM"){
    fit_iGREX <- MoM_3var(fit_init$K/G,Ka,fit_init$z,fit_init$cov_GWAS[,-1])
  }
  return(list(fit_init,fit_iGREX))
}
