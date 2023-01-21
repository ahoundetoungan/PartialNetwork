// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include "tools.h"

using namespace Rcpp;
//using namespace arma;
using namespace std;
//using namespace Eigen;
//using namespace Numer;



// update theta
void updtheta (arma::vec& theta,
               List& Vtheta,
               List& Xb,
               List& Xgamma,
               const double& sigma2,
               const double& kv,
               const double& kbeta,
               const List& Xone,
               const List& X,
               const List& Ay,
               const List& V,
               const arma::vec& invsigmathetatheta0,
               const arma::mat& invsigmatheta,
               const double& M) {
  arma::vec Amy = Ay[0];
  arma::mat Vm  = V[0];
  arma::vec VAY = Vm.t()*Amy;
  arma::mat VV = Vm.t()*Vm;
  
  for (int m(1); m<M; ++m) {
    arma::vec Amy = Ay[m];
    arma::mat Vm  = V[m];
    VAY += Vm.t()*Amy;
    VV += Vm.t()*Vm;
  }
  
  const arma::vec muhtheta = VAY + invsigmathetatheta0;
  const arma::mat sigmahtheta = arma::inv(VV + invsigmatheta);
  theta = Fmvnorm(kv,sigmahtheta*muhtheta,sigma2*sigmahtheta);
  
  for(int m(0); m<M; ++m){
    arma::mat Vm = V[m];
    arma::mat Xonem = Xone[m];
    arma::mat Xm = X[m];
    arma::vec Vmtheta = Vm*theta;
    
    // Vtheta
    Vtheta[m] = Vmtheta;
    
    //Xb
    Xb[m] = Xonem*theta.head(kbeta);
    
    //Xgamma
    Xgamma[m] = Xm*theta.tail(kv-kbeta);
  }
}

void updthetanoc (arma::vec& theta,
                  List& Vtheta,
                  const double& sigma2,
                  const double& kv,
                  const List& Ay,
                  const List& V,
                  const arma::vec& invsigmathetatheta0,
                  const arma::mat& invsigmatheta,
                  const double& M) {
  arma::vec Amy = Ay[0];
  arma::mat Vm  = V[0];
  arma::vec VAY = Vm.t()*Amy;
  arma::mat VV = Vm.t()*Vm;
  
  for (int m(1); m<M; ++m) {
    arma::vec Amy = Ay[m];
    arma::mat Vm  = V[m];
    VAY += Vm.t()*Amy;
    VV += Vm.t()*Vm;
  }
  
  const arma::vec muhtheta = VAY + invsigmathetatheta0;
  const arma::mat sigmahtheta = arma::inv(VV + invsigmatheta);
  theta = Fmvnorm(kv,sigmahtheta*muhtheta,sigma2*sigmahtheta);
  
  for(int m(0); m<M; ++m){
    arma::mat Vm = V[m];
    arma::vec Vmtheta = Vm*theta;
    
    // Vtheta
    Vtheta[m] = Vmtheta;
  }
}


// update sigma
void updsigma2 (double& sigma2,
                const arma::vec& theta,
                const double& a,
                const double& b,
                const arma::vec theta0,
                const arma::mat& invsigmatheta,
                const List& Ay,
                const List& Vtheta,
                const double& sumN,
                const double& M) {
  // ah
  double ah = 0.5*(a+sumN);
  arma::vec temp1 = theta - theta0;
  
  // bh
  double bh = b + arma::dot(temp1, invsigmatheta*temp1);
  
  arma::vec Amy = Ay[0];
  arma::vec Vmtheta  = Vtheta[0];
  arma::vec temp2 = Amy - Vmtheta;
  bh += arma::dot(temp2, temp2);
  for (int m(1); m<M; ++m) {
    arma::vec Amy = Ay[m];
    arma::vec Vmtheta  = Vtheta[m];
    temp2 = Amy - Vmtheta;
    bh += arma::dot(temp2, temp2);
  }
  bh *= 0.5;
  
  // draw sigma
  NumericVector sigma2temps = rgamma(1,ah,1/bh);
  sigma2 = 1/sigma2temps(0);
}


// update zeta
void updzeta (double& zeta,
              double& alpha,
              List& A,
              double& sumlogdetA,
              List& Ay,
              const List& Gnorm,
              const List& y,
              const double& sigma2,
              const List& Vtheta,
              const double& jumpzeta,
              double& zetaaccept,
              const double& zeta0,
              const double& invsigma2zeta,
              const Rcpp::IntegerVector N,
              const double M) {
  
  double zetastart =  R::rnorm(zeta, jumpzeta);
  double alphastart = (exp(zetastart) - 1.0)/(exp(zetastart) + 1);
  
  // Declaration of some variable 
  double logalphazeta, logalpha2zeta;
  
  //Compute logalpha2
  List Astart(M), Aystart(M);
  
  arma::mat Gm = Gnorm[0];
  int Nm = N(0);
  arma::vec Aym = Ay[0];
  arma::vec ym = y[0];
  arma::vec Vmtheta  = Vtheta[0];
  
  arma::mat Astartm = arma::eye(Nm,Nm) - alphastart*Gm;
  arma::vec Aystartm = Astartm*ym;
  Astart[0] = Astartm; // save A*
  Aystart[0] = Aystartm; // save A*y
  
  double valm, signm;
  arma::log_det(valm, signm, Astartm);
  double sumlogdetAstart = valm + log(signm); // compute sumlogdetA
  
  logalpha2zeta = arma::dot(Aystartm,Aystartm) - arma::dot(Aym,Aym) +
    2*arma::dot(Vmtheta, Aym - Aystartm);
  
  for(int m(1); m<M; ++m) {
    arma::mat Gm = Gnorm[m];
    int Nm = N(m);
    arma::vec Aym = Ay[m];
    arma::vec ym = y[m];
    arma::vec Vmtheta  = Vtheta[m];
    
    Astartm = arma::eye(Nm,Nm) - alphastart*Gm;
    Aystartm = Astartm*ym;
    Astart[m] = Astartm;
    Aystart[m] = Aystartm;
    
    double valm, signm;
    arma::log_det(valm, signm, Astartm);
    sumlogdetAstart += valm + log(signm);

    
    logalpha2zeta += arma::dot(Aystartm,Aystartm) - arma::dot(Aym,Aym) +
      2*arma::dot(Vmtheta, Aym - Aystartm);
  }
  
  logalpha2zeta *= (-.5/sigma2);
  logalpha2zeta += sumlogdetAstart - sumlogdetA + 
    (0.5*invsigma2zeta)*(pow(zeta-zeta0,2) - pow(zetastart-zeta0,2));
  
  //Compute logalpha
  logalphazeta=min(Rcpp::NumericVector::create(0, logalpha2zeta));
  
  if(unif_rand()<exp(logalphazeta)){
    zeta = zetastart;
    alpha = alphastart;
    A = Astart;
    Ay = Aystart;
    sumlogdetA = sumlogdetAstart;
    zetaaccept +=1;     //Increase acceptance number to 1 
  }
}


