// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>
#include "lmm_pxem.hpp"
#include "data_loader.hpp"
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
RcppExport SEXP lmm_pxem(const arma::vec y, const arma::mat X, const arma::mat W, const double tol, const int maxIter){

  double sigma2y, sigma2beta, loglik;
  vec beta0 =zeros<vec>(W.n_cols);
  int iter;
  mat Sigb = zeros<mat>(X.n_cols,X.n_cols);
  vec mub  = zeros<vec>(X.n_cols);

  lmm_pxem_ptr(y, W, X, tol, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);


  Rcpp::List ret;
  ret["beta0"] = beta0;
  ret["sigma2y"] = sigma2y;
  ret["sigma2beta"] = sigma2beta;
  ret["Sigb"] = Sigb;
  ret["mub"] = mub;
  ret["loglik"] = loglik;
  ret["iter"] = iter;
  return ret;
}

// [[Rcpp::export]]
RcppExport SEXP iGREX_Kg(const arma::vec& y, const arma::mat& X1, const arma::mat& X2, const arma::mat& W1, const double tol, const int maxIter) {
  int n1 = y.n_elem, n2 = X2.n_rows, p1 = X1.n_cols, p2 = X2.n_cols;
  if (p1 != p2){
    perror("The dimensions of X1 and X2 are not matched");
  }
  mat X1tmp(n1,p1);
  mat X2tmp(n2,p2);


  rowvec meanX1tmp = mean(X1, 0);
  rowvec sdX1tmp = stddev(X1, 0, 0); // see manual
  X1tmp = (X1 - repmat(meanX1tmp, n1, 1))/ repmat(sdX1tmp, n1, 1) / sqrt(p1);

  rowvec meanX2tmp = mean(X2, 0);
  rowvec sdX2tmp = stddev(X2, 0, 0); // see manual
  X2tmp = (X2 - repmat(meanX2tmp, n2, 1))/ repmat(sdX2tmp, n2, 1) / sqrt(p2);

  // step 1: LMM_PXEM
  double sigma2y, sigma2beta, loglik;
  vec beta0 =zeros<vec>(W1.n_cols);
  int iter;
  mat Sigb = zeros<mat>(p1,p1);
  vec mub  = zeros<vec>(p1);

  lmm_pxem_ptr(y, W1, X1tmp, tol, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);


  // step 2: get mat K of g-th gene by MoM
  vec X2mu = X2tmp * mub;
  mat Kg = X2mu * X2mu.t() + X2tmp * Sigb * X2tmp.t();
  mat Kg0 = X2mu * X2mu.t();
  mat weightMat = mub*mub.t() + Sigb;
  // int B = 50;
  // mat zb;
  // zb.randn(n2,B);


  Rcpp::List ret;
  ret["beta0"] = beta0;
  ret["sigma2y"] = sigma2y;
  ret["sigma2beta"] = sigma2beta;
  ret["K_g"] = Kg;
  ret["K_g0"] = Kg0;
  ret["mub"] = mub;
  ret["Sigb"] = Sigb;
  ret["loglik"] = loglik;
  ret["iter"] = iter;
  ret["weight"] = weightMat;
  return ret;
}


// [[Rcpp::export]]
RcppExport SEXP iGREX_raw(std::string prefix_eQTL_geno, std::string prefix_GWAS, std::string gene_expr, std::string cov_eQTL,
                            std::string cov_GWAS, int whCol, int bw, int subsample){ //int normalize_option = 1, int pred_option = 0){//, char* A21, char* A22){
  // normalize_option: 1. normalize each separately, 2. normalize both plink files together
  // match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
  // plink file 1: prefix_eQTL_geno; plink file 2: prefix_GWAS; expression file: gene_expr
  // covariates file for file 1: cov_eQTL; covariates file for file 2: cov_GWAS
  // pred_option :0 (no calculation for prediction) 1 (calcuation for prediction)

  List tmp = dataLoader1(prefix_eQTL_geno, prefix_GWAS, gene_expr, whCol);
  vec z = tmp["y"];
  mat expr = tmp["expr_used"];
  CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
  uvec chr_4use_r = tmp["chr_4use_r"];
  uvec bp_4use_r = tmp["bp_4use_r"];
  CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], targetID = tmp["targetID"], indiv_4use = tmp["indiv_4use"];
  vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"], ind = tmp["ind"];
  uvec idxin1 = tmp["idxin1"], idxin2 = tmp["idxin2"], idxinFile1 = tmp["idxinFile1"], idxinFile2 = tmp["idxinFile2"];

  // load the size of plink file
  string famfile1 = prefix_eQTL_geno;
  famfile1 += ".fam";
  int N1 = getLineNum(famfile1);
  string bimfile1 = prefix_eQTL_geno;
  bimfile1 += ".bim";
  long int P1 =  getLineNum(bimfile1);
  string famfile2 = prefix_GWAS;
  famfile2 += ".fam";
  int N2 = getLineNum(famfile2);
  string bimfile2 = prefix_GWAS;
  bimfile2 += ".bim";
  long int P2 =  getLineNum(bimfile2);

  long long size1 = (long long)N1 * (long long)P1;
  long long size2 = (long long)N2 * (long long)P2;

  char* X1 = new char[size1];
  char* X2 = new char[size2];
  char delimiter = '\t';

  clock_t t1 = clock();
  cout << "## Start loading genotype files 1, " ;
  readPlink2(prefix_eQTL_geno,N1, P1, X1);
  cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

  t1 = clock();
  cout << "## Start loading genotype files 2, ";
  readPlink2(prefix_GWAS,N2, P2, X2);
  cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

  cout << "## Start loading covariates files ... " << endl;
  // load covariates file w2
  mat covar2;
  if (!cov_GWAS.empty()){
    List tmp = getColNum_Header(cov_GWAS, delimiter);
    int Ncovar = tmp["columns"];
    tmp = getCovarFile(cov_GWAS, delimiter, Ncovar, N2);
    mat covartmp = tmp["covar"];
    mat w2one = ones<mat>(N2, 1);

    covar2 = join_rows(w2one, covartmp);
  }
  else {
    covar2 = ones<mat>(N2,1);
  }

  // load covariates file w1
  mat covar1;
  CharacterVector IID_w1;
  if (!cov_eQTL.empty()){
    List tmp = getColNum_Header(cov_eQTL, delimiter);

    int Ncovar = tmp["columns"];
    tmp = getCovarFile(cov_eQTL, delimiter, Ncovar, idxin1.n_elem);
    mat covartmp = tmp["covar"];
    IID_w1 = tmp["IID"];
    IntegerVector idx_tmp = match(indiv_4use, IID_w1) -1;
    uvec idx_exprcov(as<uvec>(idx_tmp));

    covar1 = join_rows(ones<mat>(idxin1.n_elem, 1), covartmp.rows(idx_exprcov));
  }
  else {
    covar1 = ones<mat>(idxin1.n_elem,1);
  }
  cout << "## End loading files ... " << endl;

  mat w1 = covar1, w2 = covar2;

  uword Ngene = lower.size();

  uvec idx;

  vec bp_4use_r_vec = conv_to<vec>::from(bp_4use_r);
  vec chr_4use_r_vec = conv_to<vec>::from(chr_4use_r);

  mat out_param = -99*ones<mat>(Ngene,5);

  uvec idx_all = zeros<uvec>(0);
  uvec idx_active_gene = zeros<uvec>(0);
  uvec g_tmp(1);
  int maxIter = 1000;
  double tol = 1e-5;
  int constr;

  mat K(N2,N2,fill::zeros);
  mat Krd(subsample,subsample,fill::zeros);
  uvec idx_subsample(subsample,fill::zeros);
  uvec sequence;

  if(subsample > 0){
    sequence = linspace<uvec>(0,N2-1,N2);
    idx_subsample = Rcpp::RcppArmadillo::sample(sequence,subsample,false);
  } else if(subsample < 0){
    perror("The value of `subsample' must be non-negative!");
  }


  // conduct fitting

  t1 = clock();
  for (uword g = 0; g < Ngene; g++){

    idx = find(bp_4use_r_vec < upper(g) + bw && bp_4use_r_vec > lower(g) - bw
                 && chr_4use_r_vec == chr_expr(g));


    out_param(g, 3) = idx.n_elem;


    if (idx.is_empty() == false){
      if (idx.n_elem > 1){
        g_tmp(0) = g;
        idx_active_gene = join_cols(idx_active_gene,g_tmp);


        if ( idx_active_gene.n_elem % 100 == 0 && idx_active_gene.n_elem != 0){
          cout << idx_active_gene.n_elem << "-th Gene starts working ..." ;
          cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
        }

        idx_all = join_cols(idx_all, idx);

        double* sub_matrix_double1 = new double[idx.n_elem * N1];
        double* sub_matrix_double2 = new double[idx.n_elem * N2];

        mat X1tmp = getSubMat(X1, N1 , P1, idxinFile1(idx)-1, sub_matrix_double1);
        X1tmp = X1tmp.rows(idxin1);
        X1tmp.replace(3, 0);
        mat X2tmp = getSubMat(X2, N2 , P2, idxinFile2(idx)-1, sub_matrix_double2);
        X2tmp.replace(3, 0);


        vec ind_idx = ind(idx);
        uvec ind_idx_1 = find(ind_idx == -1);

        X2tmp.cols(ind_idx_1) = 2 - X2tmp.cols(ind_idx_1);

        uvec idx1(1); idx1(0) = g;
        vec y = trans(expr.rows(idx1));
        y.replace(datum::nan, 9999);

        uvec idx2 = find(y != 9999);
        X1tmp = ((&X1tmp) -> rows(idx2));


        mat w1tmp = ((&w1) -> rows(idx2));

        y = y(idx2);

        rowvec meanX1tmp = mean(X1tmp, 0);
        rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
        X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1))/ repmat(sdX1tmp, X1tmp.n_rows, 1) / sqrt(X1tmp.n_cols);

        rowvec meanX2tmp = mean(X2tmp, 0);
        rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual
        X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);

        rowvec X1row = X1tmp.row(0);
        X1tmp = X1tmp.cols(find_finite(X1row));

        X2tmp = X2tmp.cols(find_finite(X1row));

        // initialize by linear mixed model for sigma2y and sigma2beta
        // use this initialization to calculate MoM estimates of h2
        double sigma2y, sigma2beta, loglik;
        vec beta0 =zeros<vec>(w1tmp.n_cols);
        int iter;
        mat Sigb = zeros<mat>(X1tmp.n_cols,X1tmp.n_cols);
        vec mub  = zeros<vec>(X1tmp.n_cols);

        // step 1: LMM_PXEM
        lmm_pxem_ptr(y, w1tmp, X1tmp, tol, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);
        out_param(g, 4) = 1/(1+sigma2y/sigma2beta); // gene-wise h2

        // step 2: get mat K of g-th gene for MoM

        vec X2mu = X2tmp * mub;
        mat Kg = X2mu * X2mu.t() + X2tmp * Sigb * X2tmp.t();
        K += Kg;

        if(subsample > 0){
          X2tmp = X2tmp.rows(idx_subsample);
          vec X2mu = X2tmp * mub;
          mat Kg = X2mu * X2mu.t() + X2tmp * Sigb * X2tmp.t();
          Krd += Kg;
        }

        //remove local variable by reset
        X1tmp.reset();
        X2tmp.reset();
        w1tmp.reset();
        Sigb.reset();
        mub.reset();
        meanX1tmp.reset();
        meanX2tmp.reset();
        sdX1tmp.reset();
        sdX2tmp.reset();
        X1row.reset();
        ind_idx.reset();
        ind_idx_1.reset();
        X2mu.reset();
        Kg.reset();
        sequence.reset();

        delete[] sub_matrix_double1;
        delete[] sub_matrix_double2;

        out_param(g, 0) = sigma2beta;
        out_param(g, 1) = sigma2y;
        out_param(g, 2) = loglik;
      }
      else{
      }
    }

    else{
    }

  }
  delete[] X1;
  delete[] X2;

  cout << "Model fitting is done in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

  cout << "Start calculting mediate Heritability" << endl;

  double sigma2g;
  double sigma2z;
  double med_H;
  double se_H;
  double se_sigma2g;
  double denom;
  mat temp;


  double trK2 = accu(pow(K,2));
  double trK = trace(K);
  double tr2K = pow(trK,2);

  mat diagN2(N2,N2,fill::eye);

  if(subsample == 0){
    denom = trK2-tr2K/N2;
    temp = K-diagN2*trK/N2;

    sigma2g = as_scalar(z.t() * temp * z) / denom;
    sigma2z = as_scalar(z.t() * (diagN2*trK2/N2-K*trK/N2) * z) / denom;

    mat Sigma = sigma2g * K + sigma2z * diagN2;
    med_H = trK * sigma2g / as_scalar(sum(pow(z,2)));
    se_sigma2g = sqrt(2*trace(z * z.t() * temp * Sigma * temp)) / denom;
    se_H = trK * se_sigma2g / as_scalar(sum(pow(z,2)));
  } else {
    double trKrd2 = accu(pow(Krd,2));
    double trKrd = trace(Krd);
    double tr2Krd = pow(trKrd,2);

    denom = trKrd2-tr2Krd/(subsample-1);
    temp = K-diagN2*trK/(N2-1);

    sigma2g = as_scalar(z.t() * temp * z) / denom * pow(subsample-1,2)/pow(N2-1,2);
    sigma2z = (as_scalar(sum(pow(z,2)))-trK*sigma2g)/(N2-1);

    mat Sigma = sigma2g * K + sigma2z * diagN2;
    med_H = trK * sigma2g / as_scalar(sum(pow(z,2)));
    se_sigma2g = sqrt(2*trace(z * z.t() * temp * Sigma * temp)) / denom;
    se_H = trK * se_sigma2g / as_scalar(sum(pow(z,2)));
  }

  // combine gene info into mat

  Rcpp::DataFrame snp_info = Rcpp::DataFrame::create(Rcpp::Named("Chr")=chr_4use_r,
                                                     Rcpp::Named("rsname")=rsname_4use_r,
                                                     Rcpp::Named("BP")=bp_4use_r);
  // Rcpp::IntegerVector(bp_4use_r_vec.begin(),bp_4use_r_vec.end())

  Rcpp::DataFrame gene_info = Rcpp::DataFrame::create(Rcpp::Named("lower")=lower,
                                                      Rcpp::Named("upper")=upper,
                                                      Rcpp::Named("genetype1")=genetype1,
                                                      Rcpp::Named("genetype2")=genetype2,
                                                      Rcpp::Named("TargetID")=targetID,
                                                      Rcpp::Named("Chr")=chr_expr);

  vec lower_active = lower(idx_active_gene);
  vec upper_active = upper(idx_active_gene);
  vec chr_active = chr_expr(idx_active_gene);

  Rcpp::DataFrame gene_info0 = Rcpp::DataFrame::create(Rcpp::Named("lower")=lower_active,
                                                      Rcpp::Named("upper")=upper_active,
                                                      Rcpp::Named("genetype1")=genetype1[as<Rcpp::IntegerVector>(wrap(idx_active_gene))],
                                                      Rcpp::Named("genetype2")=genetype2[as<Rcpp::IntegerVector>(wrap(idx_active_gene))],
                                                      Rcpp::Named("TargetID")=targetID[as<Rcpp::IntegerVector>(wrap(idx_active_gene))],
                                                      Rcpp::Named("Chr")=chr_active);

  uword Ngene_active = idx_active_gene.n_elem;
  mat out_param0 = out_param.rows(idx_active_gene);

  List out = List::create(Rcpp::Named("idx_all") = idx_all,
                          Rcpp::Named("z") = z,
                          Rcpp::Named("covar1") = w1,
                          Rcpp::Named("covar2") = w2,
                          Rcpp::Named("expr") = expr,
                          Rcpp::Named("snp_info") = snp_info,
                          Rcpp::Named("gene_info_all") = gene_info,
                          Rcpp::Named("gene_info_match") = gene_info0,
                          Rcpp::Named("out_param") = out_param0,
                          Rcpp::Named("sigma2g") = sigma2g,
                          Rcpp::Named("sigma2z") = sigma2z,
                          Rcpp::Named("se_H") = se_H,
                          Rcpp::Named("med_H") = med_H,
                          Rcpp::Named("K") = K,
                          Rcpp::Named("Krd") = Krd,
                          Rcpp::Named("subsample") = subsample,
                          Rcpp::Named("idx_subsample") = idx_subsample);

  return out;
}


// [[Rcpp::export]]
RcppExport SEXP iGREX_init(std::string prefix_eQTL_geno, std::string prefix_GWAS, std::string gene_expr, std::string cov_eQTL,
                     std::string cov_GWAS, int whCol, int bw, int subsample){ //int normalize_option = 1, int pred_option = 0){//, char* A21, char* A22){
  // normalize_option: 1. normalize each separately, 2. normalize both plink files together
  // match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
  // plink file 1: prefix_eQTL_geno; plink file 2: prefix_GWAS; expression file: gene_expr
  // covariates file for file 1: cov_eQTL; covariates file for file 2: cov_GWAS
  // pred_option :0 (no calculation for prediction) 1 (calcuation for prediction)

  List tmp = dataLoader1(prefix_eQTL_geno, prefix_GWAS, gene_expr, whCol);
  vec z = tmp["y"];
  mat expr = tmp["expr_used"];
  CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
  uvec chr_4use_r = tmp["chr_4use_r"];
  uvec bp_4use_r = tmp["bp_4use_r"];
  CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], targetID = tmp["targetID"], indiv_4use = tmp["indiv_4use"];
  vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"], ind = tmp["ind"];
  uvec idxin1 = tmp["idxin1"], idxin2 = tmp["idxin2"], idxinFile1 = tmp["idxinFile1"], idxinFile2 = tmp["idxinFile2"];

  // load the size of plink file
  string famfile1 = prefix_eQTL_geno;
  famfile1 += ".fam";
  int N1 = getLineNum(famfile1);
  string bimfile1 = prefix_eQTL_geno;
  bimfile1 += ".bim";
  long int P1 =  getLineNum(bimfile1);
  string famfile2 = prefix_GWAS;
  famfile2 += ".fam";
  int N2 = getLineNum(famfile2);
  string bimfile2 = prefix_GWAS;
  bimfile2 += ".bim";
  long int P2 =  getLineNum(bimfile2);

  long long size1 = (long long)N1 * (long long)P1;
  long long size2 = (long long)N2 * (long long)P2;

  char* X1 = new char[size1];
  char* X2 = new char[size2];
  char delimiter = '\t';

  clock_t t1 = clock();
  cout << "## Start loading genotype files 1, " ;
  readPlink2(prefix_eQTL_geno,N1, P1, X1);
  cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

  t1 = clock();
  cout << "## Start loading genotype files 2, ";
  readPlink2(prefix_GWAS,N2, P2, X2);
  cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

  cout << "## Start loading covariates files ... " << endl;
  // load covariates file w2

  mat covar2;
  if (!cov_GWAS.empty()){
    List tmp = getColNum_Header(cov_GWAS, delimiter);
    int Ncovar = tmp["columns"];
    tmp = getCovarFile(cov_GWAS, delimiter, Ncovar, N2);
    mat covartmp = tmp["covar"];
    mat w2one = ones<mat>(N2, 1);

    covar2 = join_rows(w2one, covartmp);
  }
  else {
    covar2 = ones<mat>(N2,1);
  }

  // load covariates file w1
  mat covar1;
  CharacterVector IID_w1;
  if (!cov_eQTL.empty()){

    List tmp = getColNum_Header(cov_eQTL, delimiter);

    int Ncovar = tmp["columns"];
    tmp = getCovarFile(cov_eQTL, delimiter, Ncovar, idxin1.n_elem);
    mat covartmp = tmp["covar"];
    IID_w1 = tmp["IID"];
    IntegerVector idx_tmp = match(indiv_4use, IID_w1) -1;
    uvec idx_exprcov(as<uvec>(idx_tmp));

    covar1 = join_rows(ones<mat>(idxin1.n_elem, 1), covartmp.rows(idx_exprcov));
  }
  else {
    covar1 = ones<mat>(idxin1.n_elem,1);
  }
  cout << "Number of covariates for eQTL: " << covar1.n_cols << endl;
  cout << "Number of covariates for GWAS: " << covar2.n_cols << endl;
  cout << "## End loading files ... " << endl;

  mat w1 = covar1, w2 = covar2;

  uword Ngene = lower.size();

  uvec idx;

  vec bp_4use_r_vec = conv_to<vec>::from(bp_4use_r);
  vec chr_4use_r_vec = conv_to<vec>::from(chr_4use_r);

  mat out_param = -99*ones<mat>(Ngene,5);

  uvec idx_all = zeros<uvec>(0);
  uvec idx_active_gene = zeros<uvec>(0);
  uvec g_tmp(1);
  int maxIter = 1000;
  double tol = 1e-5;
  int constr;

  mat K(N2,N2,fill::zeros);
  mat K0(N2,N2,fill::zeros);
  uvec idx_subsample(subsample,fill::zeros);
  uvec sequence;

  if(subsample > 0){
    sequence = linspace<uvec>(0,N2-1,N2);
    idx_subsample = Rcpp::RcppArmadillo::sample(sequence,subsample,false);
  } else if(subsample < 0){
    perror("The value of `subsample' must be non-negative!");
  }


  // conduct fitting

  t1 = clock();
  for (uword g = 0; g < Ngene; g++){

    idx = find(bp_4use_r_vec < upper(g) + bw && bp_4use_r_vec > lower(g) - bw
                 && chr_4use_r_vec == chr_expr(g));


    out_param(g, 3) = idx.n_elem;


    if (idx.is_empty() == false){
      if (idx.n_elem > 1){
        g_tmp(0) = g;
        idx_active_gene = join_cols(idx_active_gene,g_tmp);


        if ( idx_active_gene.n_elem % 100 == 0 && idx_active_gene.n_elem != 0){
          cout << idx_active_gene.n_elem << "-th Gene starts working ..." ;
          cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
        }

        // idx_all = join_cols(idx_all, idx);

        double* sub_matrix_double1 = new double[idx.n_elem * N1];
        double* sub_matrix_double2 = new double[idx.n_elem * N2];

        mat X1tmp = getSubMat(X1, N1 , P1, idxinFile1(idx)-1, sub_matrix_double1);
        X1tmp = X1tmp.rows(idxin1);
        X1tmp.replace(3, 0);
        mat X2tmp = getSubMat(X2, N2 , P2, idxinFile2(idx)-1, sub_matrix_double2);
        X2tmp.replace(3, 0);


        vec ind_idx = ind(idx);
        uvec ind_idx_1 = find(ind_idx == -1);

        X2tmp.cols(ind_idx_1) = 2 - X2tmp.cols(ind_idx_1);

        uvec idx1(1); idx1(0) = g;
        vec y = trans(expr.rows(idx1));
        y.replace(datum::nan, 9999);

        uvec idx2 = find(y != 9999);
        X1tmp = ((&X1tmp) -> rows(idx2));


        mat w1tmp = ((&w1) -> rows(idx2));

        y = y(idx2);

        rowvec meanX1tmp = mean(X1tmp, 0);
        rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
        X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1))/ repmat(sdX1tmp, X1tmp.n_rows, 1) / sqrt(X1tmp.n_cols);

        rowvec meanX2tmp = mean(X2tmp, 0);
        rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual
        X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);

        rowvec X1row = X1tmp.row(0);
        X1tmp = X1tmp.cols(find_finite(X1row));

        X2tmp = X2tmp.cols(find_finite(X1row));

        idx = idx.elem(find_finite(X1row));
        idx_all = join_cols(idx_all, idx);

        // initialize by linear mixed model for sigma2y and sigma2beta
        // use this initialization to calculate MoM estimates of h2
        double sigma2y, sigma2beta, loglik;
        vec beta0 =zeros<vec>(w1tmp.n_cols);
        int iter;
        mat Sigb = zeros<mat>(X1tmp.n_cols,X1tmp.n_cols);
        vec mub  = zeros<vec>(X1tmp.n_cols);

        // step 1: LMM_PXEM
        lmm_pxem_ptr(y, w1tmp, X1tmp, tol, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);
        out_param(g, 4) = 1/(1+sigma2y/sigma2beta); // gene-wise h2

        // step 2: get mat K of g-th gene for MoM
        if(subsample>0){
          X2tmp = X2tmp.rows(idx_subsample);
        }
        vec X2mu = X2tmp * mub;
        mat Kg = X2mu * X2mu.t() + X2tmp * Sigb * X2tmp.t();
        K += Kg;

        Kg = X2mu * X2mu.t();
        K0 += Kg;

        //remove local variable by reset
        X1tmp.reset();
        X2tmp.reset();
        w1tmp.reset();
        Sigb.reset();
        mub.reset();
        meanX1tmp.reset();
        meanX2tmp.reset();
        sdX1tmp.reset();
        sdX2tmp.reset();
        X1row.reset();
        ind_idx.reset();
        ind_idx_1.reset();
        X2mu.reset();
        Kg.reset();
        sequence.reset();

        delete[] sub_matrix_double1;
        delete[] sub_matrix_double2;

        out_param(g, 0) = sigma2beta;
        out_param(g, 1) = sigma2y;
        out_param(g, 2) = loglik;
      }
      else{
      }
    }

    else{
    }

  }
  delete[] X1;
  delete[] X2;

  cout << "Model fitting is done in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

  // combine gene info into mat

  Rcpp::DataFrame snp_info = Rcpp::DataFrame::create(Rcpp::Named("Chr")=chr_4use_r,
                                                     Rcpp::Named("rsname")=rsname_4use_r,
                                                     Rcpp::Named("BP")=bp_4use_r);
  // Rcpp::IntegerVector(bp_4use_r_vec.begin(),bp_4use_r_vec.end())

  Rcpp::DataFrame gene_info = Rcpp::DataFrame::create(Rcpp::Named("lower")=lower,
                                                      Rcpp::Named("upper")=upper,
                                                      Rcpp::Named("genetype1")=genetype1,
                                                      Rcpp::Named("genetype2")=genetype2,
                                                      Rcpp::Named("TargetID")=targetID,
                                                      Rcpp::Named("Chr")=chr_expr);

  vec lower_active = lower(idx_active_gene);
  vec upper_active = upper(idx_active_gene);
  vec chr_active = chr_expr(idx_active_gene);

  Rcpp::DataFrame gene_info0 = Rcpp::DataFrame::create(Rcpp::Named("lower")=lower_active,
                                                       Rcpp::Named("upper")=upper_active,
                                                       Rcpp::Named("genetype1")=genetype1[as<Rcpp::IntegerVector>(wrap(idx_active_gene))],
                                                                                         Rcpp::Named("genetype2")=genetype2[as<Rcpp::IntegerVector>(wrap(idx_active_gene))],
                                                                                                                           Rcpp::Named("TargetID")=targetID[as<Rcpp::IntegerVector>(wrap(idx_active_gene))],
                                                                                                                                                           Rcpp::Named("Chr")=chr_active);
  Rcpp::DataFrame idxinFile = Rcpp::DataFrame::create(Rcpp::Named("idxinFile1")=idxinFile1,
                                                      Rcpp::Named("idxinFile2")=idxinFile2);

  uword Ngene_active = idx_active_gene.n_elem;
  mat out_param0 = out_param.rows(idx_active_gene);

  List out = List::create(Rcpp::Named("idx_all") = idx_all,
                          Rcpp::Named("idxinFile") = idxinFile,
                          Rcpp::Named("ind") = ind,
                          Rcpp::Named("z") = z,
                          Rcpp::Named("covar1") = w1,
                          Rcpp::Named("covar2") = w2,
                          Rcpp::Named("expr") = expr,
                          Rcpp::Named("snp_info") = snp_info,
                          Rcpp::Named("gene_info_all") = gene_info,
                          Rcpp::Named("gene_info_match") = gene_info0,
                          Rcpp::Named("out_param") = out_param0,
                          Rcpp::Named("K") = K,
                          Rcpp::Named("K0") = K0,
                          Rcpp::Named("subsample") = subsample,
                          Rcpp::Named("idx_subsample") = idx_subsample);

  return out;
}


/*** R
timesTwo(42)
*/
