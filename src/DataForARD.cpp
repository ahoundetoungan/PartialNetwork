// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define NDEBUG 1

using namespace Rcpp;
using namespace arma;
using namespace std;

//The fonction Prob computes the link probabilities of vertices
// [[Rcpp::export]]
arma::mat Prob(arma::vec& nu, arma::vec& d, double& zeta, arma::mat& z){
  int N=z.n_rows;
  // compute Probabilities
  arma::mat prob = zeta*z*z.t();
  
  prob += arma::repmat(nu,1,N) + arma::repmat(nu.t(),N,1);
  prob = arma::exp(prob);
  prob.diag()=zeros(N);   //zero on the diagonal
  prob*=(arma::sum(d)/arma::accu(prob));
  
  arma::umat tempbin = (prob<1);
  arma::mat tempone = arma::mat(N,N,fill::ones);
  prob = tempbin%prob + (1-tempbin)%tempone;
  
  return prob;
}




//The function Graph conputes the Adjacency Matrix
// [[Rcpp::export]]
arma::umat Graph(arma::mat& prob){
  int N=prob.n_rows;
  arma::mat matunif(N,N,fill::randu);
  //mat G(N,N);
  
 // for(int i(0);i<N;++i){
  //  for(int j(0);j<N;++j){
  //    G(i,j)=rbinom(1,1,prob(i,j))(0);
  //  }
  //}
  return (matunif<prob);
  //return G;
}