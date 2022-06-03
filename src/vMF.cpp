// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define NDEBUG 1

using namespace Rcpp;
using namespace arma;
using namespace std;


/*
 *  Estimate graph parameters Mccormick and Zheng (2015)
 *  
 *  To sample from von Mises Fisher Distribution
 *  We follow Wood, A. T. (1994). Simulation of the von Mises Fisher distribution.
 *  Communications in statistics-simulation and computation, 23(1), 157-164.
 *  
 *  In a first step we built a function rw to compute W; that is, 
 *  perform steps 0, 1 et 2.
 *  
 *  The function needs the sample size "size", lambda as intensity parameter,
 *  d as dimension p-1 et a Vector W (void) of dimension "size" to save the
 *  values corresponding to each sampling.
 *  
 *  In the Metropolis steps, one zi is draw for each i and for each iteration
 *  Thus we build a fast version of previous functions for size = 1.
 *  That allows to avoid the loop over n and is faster.
 *  This function is not exported because it is used only in the C++ code.
 */

// rw is a machin function used to sample from vMF distribution
void rw(const int& size, const double& lambda, const int& d, arma::vec& W){
  // Step 0
  // Algebraically equivalent to
  // (-2. * l + sqrt(4. * l * l + d * d)) / d
  // but numerically more stable. See 
  // Hornik, K., & Grün, B. (2014). movMF: An R package for fitting mixtures
  // of von Mises-Fisher distributions. Journal of Statistical Software, 
  // 58(10), 1-31.
  double b = d/ (sqrt(4. * lambda * lambda + d*d) + 2. * lambda);
  double x = (1. - b) / (1. + b);
  double c = lambda * x + d * log(1. - x * x);
  
  // Step 1
  // Let's declare the variables we will use
  double w, Z, U;   // distinguish w from W. W is a vector of w
  
  // Start the loop
  for(int i(0);i<size;++i){
    Step1: Z = (rbeta(1,d/2.,d/2.))(0);
    w = (1.-(1.+b)*Z)/(1.-(1.-b)*Z);
    U = (runif(1,0.,1.))(0);
    
    //step 2
    if(lambda*w+d*log(1-x*w)-c < log(U)){goto Step1;}
    W(i)=w;
  }
}


// The following function complete the algorithm by performing the step 4 and
// the rotation toward the mean directional.
// It needs the sample size and theta as intensity parameter x mean directional

// [[Rcpp::export]]
arma::mat rvMFcpp(const int& size,const arma::vec& theta){
  int p=theta.n_rows;            //hypersphere dimension
  double lambda=norm(theta);     //intensity parameter
  arma::mat X;                         //Output matrix
  
  // if lambda=0 sample uniform; that is normal/norm
  if(lambda==0){
    NumericMatrix Xtemp(size,p,rnorm(size*p).begin());
    X=normalise(as<arma::mat>(Xtemp),2,1);  //nomalize rows by their norm   
  }
  else{
    double d=p-1;
    // Compute W
    arma::vec W(size);        //Void W
    rw(size, lambda, d, W);   //Fill W using rw
    //arma::mat Wplus=repmat(W,1,d);  //Reshape to [W W W ... W] of dimension (n,d) * remplaced using each_col
    //mean direction parameter
    arma::vec mu=theta/lambda;           
    // Necessary variables declaration
    NumericMatrix Vtemp(size,d,rnorm(size*d).begin());
    arma::mat V  = normalise(as<arma::mat>(Vtemp),2,1);
    arma::mat X1 = V.each_col()%sqrt(1 - pow(W, 2));
    X=join_rows(X1,W);
    //Rotation
    // To get samples from a vMF distribution with arbitrary mean direction
    // parameter µ, X is multiplied from the right with a matrix where the
    // first (m − 1) columns consist of unitary basis vectors of
    // the subspace orthogonal to µ and the last column is equal to µ. See
    // Hornik, K., & Grün, B. (2014). movMF: An R package for fitting mixtures
    // of von Mises-Fisher distributions. Journal of Statistical Software, 
    // 58(10), 1-31.
    arma::mat Q,R;
    qr(Q, R, mu);    //QR decomposition to get subsâce orthogonal to µ
    IntegerVector seqcol=seq_len(d);
    Q=Q.cols(as<arma::uvec>(seqcol));
    Q=join_rows(Q,mu);
    X=X*Q.t();
  }
  return X;
}

void rwone(const double& lambda, const int& d, double& w){
  // Step 0
  // Algebraically equivalent to
  // (-2. * l + sqrt(4. * l * l + d * d)) / d
  // in the reference, but numerically more stable:
  double b = d/ (sqrt(4. * lambda * lambda + d*d) + 2. * lambda);
  double x = (1. - b) / (1. + b);
  double c = lambda * x + d * log(1. - x * x);
  
  // Step 1
  // Let's declare the variables we will use
  double Z, U;   // distinguish w from W. W is a vector of w
  
  // Start the loop
  Step1: Z = (rbeta(1,d/2.,d/2.))(0);
  w = (1.-(1.+b)*Z)/(1.-(1.-b)*Z);
  U = (runif(1,0.,1.))(0);
  // step 2
  if(lambda*w+d*log(1-x*w)-c < log(U)){goto Step1;}
}

arma::mat rvMFone(const int p, const arma::vec& theta){  //p is the hypersphere dimension
  double lambda=norm(theta);  // intensity parameter
  arma::mat X(1,p);                 // The ouptup matrix
  // if lambda=0 sample uniform; that is normal/norm
  if(lambda==0){
    X=(normalise(as<arma::vec>(rnorm(p,0,1)))).t();
  }
  else{
    double d=p-1;
    // compute w 
    double w;
    rwone(lambda, d, w);
    //mean direction parameter
    arma::vec mu=theta/lambda;           
    // Necessary variables declaration
    arma::vec V=normalise(as<arma::vec>(rnorm(d,0,1)));
    arma::vec X1=sqrt(1-w*w)*V; 
    arma::vec X2(1); X2(0)=w;
    X=join_rows(X1.t(),X2);
    arma::mat Q,R;
    qr(Q, R, mu);    //QR decomposition to get subsâce orthogonal to µ
    IntegerVector seqcol=seq_len(d);
    Q=Q.cols(as<arma::uvec>(seqcol));
    Q=join_rows(Q,mu);
    X=X*Q.t();
  }
  return X;
}


// log of the normalization constant as function cpvMF
// It needs the hyperparameter dimension p and the intensity parameters k
// [[Rcpp::export]]
double logCpvMFcpp(const int& p, const double& k){
  if(k==0){ /*If k=0 return 1*/  return 0;}
  return (p/2.0 - 1.0)*log(k/2.0) - lgammal(p/2.0) - log(R::bessel_i((double) k, p/2.0 - 1, 2)) - k;
}

// Even if this is not used in this program, we want also the package provides
// a way to compute the von Mises-Fisher density.
// The function needs a matrix of points z for at which the density will be
// computed and theta = intensity paramter x mean directional

// [[Rcpp::export]]
NumericVector dvMFcpp(const arma::mat& z, const arma::vec& theta, const bool& logp = false){
  NumericVector logdens = wrap(logCpvMFcpp(z.n_cols, norm(theta)) + z*theta);
  if (logp) {
    return logdens;
  }
  return exp(logdens);
}

