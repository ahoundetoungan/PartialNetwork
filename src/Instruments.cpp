// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//Compute instruments
//If y is entered the function computes a list of 2
//[1] = a list such that the entry s is the sth draw of WX where
//WX is a cube where WX[,,p] = G^p%*%X
//[2] = a matrix G%*%y
//If y is not entered the function computes only the entry [1]
//By default S=2 et p=1

// [[Rcpp::export]]
List instruments(const arma::mat& distr, arma::mat& X, Nullable<arma::vec&> y=R_NilValue, const unsigned int& S=2, const unsigned int& pow = 1){
  const int N = X.n_rows;
  const int Kx = X.n_cols;
  
  if(y.isNotNull()){
    List output = List::create(Named("GX"),Named("GY"));
    
    mat matunif(N,N,fill::randu);
    arma::mat Gt = arma::normalise(conv_to<mat>::from((matunif < distr)),1,1);
    arma::vec Y = as<arma::vec>(y);
    output(1) = Gt*Y;
    
    List LGX(S);
    arma::cube GX(N,Kx,pow);
    GX.slice(0) =  Gt*X;
    for(int p(1); p<pow; ++p){
      GX.slice(p) =  Gt*GX.slice(p-1);
    }
    LGX(0) = GX;
    
    for(int s(1); s<S; ++s){
      mat matunif(N,N,fill::randu);
      Gt = arma::normalise(conv_to<mat>::from((matunif < distr)),1,1);

      GX.slice(0) =  Gt*X;
      for(int p(1); p<pow; ++p){
        GX.slice(p) =  Gt*GX.slice(p-1);
      }
      LGX(s) = GX;
    }
    
    output(0) = LGX;
    
    return output;
  }
  else{
    List output = List::create(Named("GX"));
    
    mat matunif(N,N,fill::randu);
    arma::mat Gt = arma::normalise(conv_to<mat>::from((matunif < distr)),1,1);
    
    List LGX(S);
    arma::cube GX(N,Kx,pow);
    GX.slice(0) =  Gt*X;
    for(int p(1); p<pow; ++p){
      GX.slice(p) =  Gt*GX.slice(p-1);
    }
    LGX(0) = GX;
    
    for(int s(1); s<S; ++s){
      mat matunif(N,N,fill::randu);
      Gt = arma::normalise(conv_to<mat>::from((matunif < distr)),1,1);
      
      GX.slice(0) =  Gt*X;
      for(int p(1); p<pow; ++p){
        GX.slice(p) =  Gt*GX.slice(p-1);
      }
      LGX(s) = GX;
    }
    
    output(0) = LGX;
    
    return output;
  }
}
