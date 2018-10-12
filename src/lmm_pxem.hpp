#ifndef lmm_pxem_hpp
#define lmm_pxem_hpp

#include <RcppArmadillo.h>
//#include <armadillo>
//#include <Rcpp.h>
//#include <omp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

void lmm_pxem_ptr(const arma::vec& y, const arma::mat& w, const arma::mat& x, const double& tol, const int& maxIter,
              double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max,
              int& iteration, arma::mat& Sigb, arma::vec& mub);

#endif /* lmm_pxem_hpp */
