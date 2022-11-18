#include <RcppArmadillo.h>
#define NDEBUG

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
double f_rcpoisson(const arma::vec& beta, 
                   const arma::vec& y,
                   const arma::mat& X,
                   const Rcpp::LogicalVector& censure,
                   const arma::vec& lcensure,
                   const arma::vec& rcensure){
  // cout<<"beta: "<<beta.t()<<endl;
  arma::vec eXb = exp(X*beta);
  int n         = y.n_elem;
  double llh = 0;
  for(int i(0); i < n; ++ i){
    if(censure(i)){
      double l2 = R::ppois(rcensure(i), eXb(i), true, true);
      double l1 = R::ppois(lcensure(i), eXb(i), true, true);
      llh      += (l1 + log(exp(l2 - l1) - 1));
    } else {
      llh      += R::dpois(y(i), eXb(i), true);
    }
  }
  return -llh;
}

//[[Rcpp::export]]
Rcpp::NumericVector expy(const arma::vec& beta, 
                         const arma::mat& X,
                         const arma::vec& lcensure,
                         const arma::vec& rcensure){
  // cout<<"beta: "<<beta.t()<<endl;
  arma::vec Xb    = X*beta;
  arma::vec eXb   = exp(Xb);
  int n           = X.n_rows;
  Rcpp::NumericVector out(n);
  for(int i(0); i < n; ++ i){
    double lpa    = R::ppois(lcensure(i), eXb(i), true, true);
    double lpamo  = R::ppois(lcensure(i) - 1, eXb(i), true, true);
    double lpb    = R::ppois(rcensure(i), eXb(i), true, true);
    double lpbmo  = R::ppois(rcensure(i) - 1, eXb(i), true, true);
    out(i)        = Xb(i) + lpamo - lpa + log(exp(lpbmo - lpamo) - 1) - log(exp(lpb - lpa) - 1);
  }
  out.attr("dim") = R_NilValue;
  return exp(out);
}