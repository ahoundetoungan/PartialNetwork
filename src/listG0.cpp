// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
List flistGnorm1 (List& dnetwork,
                  arma::vec y, 
                  arma::mat Xone,
                  arma::mat X,
            const int& M) {
  arma::vec N(M);
  List Gnorm(M), ly(M), lXone(M), lX(M);
  int r2          = -1;
  int r1;
  for (int m(0); m < M; ++ m) {
    arma::mat dnm = dnetwork(m);
    int Nm        = dnm.n_rows;
    r2           += Nm;
    r1            = r2 - Nm + 1;
    arma::mat matunif(Nm, Nm, arma::fill::randu);
    arma::mat Gnm = arma::normalise(conv_to<mat>::from((matunif < dnm)),1,1);
    N(m)          = Nm;
    Gnorm(m)      = Gnm;
    ly(m)         = y.subvec(r1, r2);
    lXone(m)      = Xone.rows(r1, r2);
    lX(m)         = X.rows(r1, r2);
  }
  
  return List::create(Named("N")     = N, 
                      Named("G")     = Gnorm, 
                      Named("ly")    = ly,
                      Named("lXone") = lXone,
                      Named("lX")    = lX);
}


//[[Rcpp::export]]
List flistGnorm2 (List& dnetwork,
                  List& Gnorm,
                  arma::vec y, 
                  arma::mat Xone,
                  arma::mat X,
             const int& M) {
  arma::vec N(M);
  List ly(M), lXone(M), lX(M);
  int r2          = -1;
  int r1;
  for (int m(0); m < M; ++ m) {
    arma::mat dnm = dnetwork(m);
    arma::mat Gnm = Gnorm(m);
    int Nm        = dnm.n_rows;
    r2           += Nm;
    r1            = r2 - Nm + 1;
    Gnm           = arma::normalise(arma::ceil(Gnm),1,1);
    N(m)          = Nm;
    Gnorm(m)      = Gnm;
    ly(m)         = y.subvec(r1, r2);
    lXone(m)      = Xone.rows(r1, r2);
    lX(m)         = X.rows(r1, r2);
  }
  
  return List::create(Named("N")     = N, 
                      Named("G")     = Gnorm, 
                      Named("ly")    = ly,
                      Named("lXone") = lXone,
                      Named("lX")    = lX);
}


//[[Rcpp::export]]
List flistGnorm1nc (List& dnetwork,
                   arma::vec y, 
                   arma::mat Xone,
                   const int& M) {
  arma::vec N(M);
  List Gnorm(M), ly(M), lXone(M), lX(M);
  int r2          = -1;
  int r1;
  for (int m(0); m < M; ++ m) {
    arma::mat dnm = dnetwork(m);
    int Nm        = dnm.n_rows;
    r2           += Nm;
    r1            = r2 - Nm + 1;
    arma::mat matunif(Nm, Nm, arma::fill::randu);
    arma::mat Gnm = arma::normalise(conv_to<mat>::from((matunif < dnm)),1,1);
    N(m)          = Nm;
    Gnorm(m)      = Gnm;
    ly(m)         = y.subvec(r1, r2);
    lXone(m)      = Xone.rows(r1, r2);
  }
  
  return List::create(Named("N")     = N, 
                      Named("G")     = Gnorm, 
                      Named("ly")    = ly,
                      Named("lXone") = lXone);
}


//[[Rcpp::export]]
List flistGnorm2nc (List& dnetwork,
                  List& Gnorm,
                  arma::vec y, 
                  arma::mat Xone,
                  const int& M) {
  arma::vec N(M);
  List ly(M), lXone(M), lX(M);
  int r2          = -1;
  int r1;
  for (int m(0); m < M; ++ m) {
    arma::mat dnm = dnetwork(m);
    arma::mat Gnm = Gnorm(m);
    int Nm        = dnm.n_rows;
    r2           += Nm;
    r1            = r2 - Nm + 1;
    Gnm           = arma::normalise(arma::ceil(Gnm),1,1);
    N(m)          = Nm;
    Gnorm(m)      = Gnm;
    ly(m)         = y.subvec(r1, r2);
    lXone(m)      = Xone.rows(r1, r2);
  }
  
  return List::create(Named("N")     = N, 
                      Named("G")     = Gnorm, 
                      Named("ly")    = ly,
                      Named("lXone") = lXone);
}
