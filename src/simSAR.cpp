/* GENERAL NOTATIONS
 * y       : is the vector of the outcome values
 * X       : is the matrix of the explanatory variables. Add the intercept 
 *           if it is included in the model. The intercept will not be added automatically. 
 * G       : is the network matrix List. That is G[r] is the subnetwork of the group r. 
 *           Gs(i,j) = measures the intensity of the outgoing link from i to j. 
 * igroup  : is the matrix of groups indexes. The data should be organized by group. 
 *           igroup[s,] is a 2-dimension vector ans gives the first and the last rows
 *           for the group s.
 * ngroup  : is the number of groups.
 * n       : The sample size.
 * nvec    : Vector of agents in each subnet
 */

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD
#endif
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

//[[Rcpp::export]]
Rcpp::List fylim(const std::vector<Eigen::MatrixXd>& G,
                 const Eigen::VectorXd& talpha,
                 const Eigen::ArrayXi& igroup,
                 const Eigen::ArrayXi& nvec,
                 const int& ngroup,
                 const double& lambda,
                 const int& n,
                 const int& nthreads) {
  Eigen::VectorXd y(n), Gy(n);
  //loop over group
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for (int m = 0; m < ngroup; ++ m) {
    int nm(nvec(m));
    Eigen::MatrixXd Am(-lambda * G[m]);
    Am.diagonal().array() += 1;
    y.segment(igroup(m), nm)  = Am.colPivHouseholderQr().solve(talpha.segment(igroup(m), nm));
    Gy.segment(igroup(m), nm) = G[m] * y.segment(igroup(m), nm);
  }
#else
  for (int m = 0; m < ngroup; ++ m) {
    int nm(nvec(m));
    Eigen::MatrixXd Am(-lambda * G[m]);
    Am.diagonal().array() += 1;
    y.segment(igroup(m), nm)  = Am.colPivHouseholderQr().solve(talpha.segment(igroup(m), nm));
    Gy.segment(igroup(m), nm) = G[m] * y.segment(igroup(m), nm);
  }
#endif
  return Rcpp::List::create(Rcpp::_["y"] = y, Rcpp::_["Gy"] = Gy);
}
