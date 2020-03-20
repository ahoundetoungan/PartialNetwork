// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//Compute instruments
//retuns list a S: number of replication
//Each component of the list has 3 element if y is given or 2 otherwise
//[1] = G1y  with first draw
//[2] = W1X is a cube where WX[,,p] = G^p%*%X  with first draw
//[2] = W2X  with second draw

// [[Rcpp::export]]
List instruments1(const arma::mat& dnetwork,
                  arma::mat& X,
                  arma::vec& y,
                  const int& S,
                  const int& pow){
  const int N        = X.n_rows;
  const int Kx       = X.n_cols;
  List output(S);
  
  
  arma::mat Gt1, Gt2;
  arma::vec G1Y;
  NumericVector G1y;
  arma::cube G1X(N,Kx,pow);
  arma::cube G2X(N,Kx,pow);
  List tmp         =  List::create(Named("G1y"), Named("G1X"), Named("G2X"));
  for(int s(0); s<S; ++s){
    // first draw
    mat matunif1(N,N,fill::randu);
    Gt1            = arma::normalise(conv_to<mat>::from((matunif1 < dnetwork)),1,1);
    
    G1Y            = Gt1*y;   
    G1y            = wrap(G1Y); G1y.attr("dim") = R_NilValue;
    tmp(0)         = G1y;
    G1X.slice(0)   =  Gt1*X;
    for(int p(1); p<pow; ++p){
      G1X.slice(p) =  Gt1*G1X.slice(p-1);
    }
    tmp(1)         = G1X;
    
    // second draw
    mat matunif2(N,N,fill::randu);
    Gt2 = arma::normalise(conv_to<mat>::from((matunif2 < dnetwork)),1,1);
    
    G2X.slice(0)   =  Gt2*X;
    for(int p(1); p<pow; ++p){
      G2X.slice(p) =  Gt2*G2X.slice(p-1);
    }
    
    tmp(2)         = G2X;
    output(s)      = tmp;
  }
  
  return output;
}



// [[Rcpp::export]]
List instruments2(const arma::mat& dnetwork,
                  arma::mat& X,
                  const int& S,
                  const int& pow){
  const int N        = X.n_rows;
  const int Kx       = X.n_cols;
  List output(S);
  
  
  arma::mat Gt1, Gt2;
  arma::cube G1X(N,Kx,pow);
  arma::cube G2X(N,Kx,pow);
  List tmp         =  List::create(Named("G1X"), Named("G2X"));
  for(int s(0); s<S; ++s){
    // first draw
    mat matunif1(N,N,fill::randu);
    Gt1            = arma::normalise(conv_to<mat>::from((matunif1 < dnetwork)),1,1);
    
    G1X.slice(0)   =  Gt1*X;
    for(int p(1); p<pow; ++p){
      G1X.slice(p) =  Gt1*G1X.slice(p-1);
    }
    tmp(0)         = G1X;
    
    // second draw
    mat matunif2(N,N,fill::randu);
    Gt2 = arma::normalise(conv_to<mat>::from((matunif2 < dnetwork)),1,1);
    
    G2X.slice(0)   =  Gt2*X;
    for(int p(1); p<pow; ++p){
      G2X.slice(p) =  Gt2*G2X.slice(p-1);
    }
    
    tmp(1)         = G2X;
    output(s)      = tmp;
  }
  
  
  return output;
}
