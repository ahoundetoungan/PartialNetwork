// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define NDEBUG 1

using namespace Rcpp;
using namespace arma;
using namespace std;
// functions in this file compule explanatory variables from the network distribution and
// the init value of G in necessary. 

/////////// WITH CONTEXTUAL EFFECTS
////G will be computed because was not given
//[[Rcpp::export]]
List flistGnorm1 (List& dnetwork,
                  arma::vec& y, 
                  arma::mat& Xone,
                  arma::mat& X,
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


////G will not be computed because was given
//[[Rcpp::export]]
List flistGnorm2 (List& G,
                  arma::vec& y, 
                  arma::mat& Xone,
                  arma::mat& X,
                  const int& M) {
  arma::vec N(M);
  List Gnorm(M), ly(M), lXone(M), lX(M);
  int r2          = -1;
  int r1;
  for (int m(0); m < M; ++ m) {
    arma::mat Gnm = G(m);
    int Nm        = Gnm.n_rows;
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

/////////// WITH CONTEXTUAL EFFECTS
////G will be computed because was not given
//[[Rcpp::export]]
List flistGnorm1nc (List& dnetwork,
                    arma::vec& y, 
                    arma::mat& Xone,
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

////G will not be computed because was given
//[[Rcpp::export]]
List flistGnorm2nc (List& G,
                    arma::vec& y, 
                    arma::mat& Xone,
                    const int& M) {
  arma::vec N(M);
  List Gnorm(M), ly(M), lXone(M), lX(M);
  int r2          = -1;
  int r1;
  for (int m(0); m < M; ++ m) {
    arma::mat Gnm = G(m);
    int Nm        = Gnm.n_rows;
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


// This function simulate a G from a given dnetwork

/////////// WITH CONTEXTUAL EFFECTS
////G will be computed because was not given
//[[Rcpp::export]]
List simG (List& dnetwork,
           const arma::vec& N,
           const int& M) {
  List out(M);
    for (int m(0); m < M; ++ m) {
      int Nm        = N(m);
      arma::mat dnm = dnetwork(m);
      arma::mat matunif(Nm, Nm, arma::fill::randu);
      //out[m]        = arma::normalise(conv_to<mat>::from((matunif < dnm)),1,1);
      dnm           = arma::conv_to<mat>::from((matunif < dnm));
      
      dnm.diag()    = arma::zeros(Nm);
      out[m]        = dnm;
    }
    
    return out;
}


//[[Rcpp::export]]
List simGnorm (List& dnetwork,
              const arma::vec& N,
              const int& M) {
  List out(M);
    for (int m(0); m < M; ++ m) {
      int Nm        = N(m);
      arma::mat dnm = dnetwork(m);
      arma::mat matunif(Nm, Nm, arma::fill::randu);
      dnm           = arma::normalise(arma::conv_to<mat>::from((matunif < dnm)),1,1);
      //out[m]        = arma::conv_to<mat>::from((matunif < dnm));
      
      dnm.diag()    = arma::zeros(Nm);
      out[m]        = dnm;
    }
    
    return out;
}
