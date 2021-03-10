// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppProgress)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppEigen.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "update.h"
#include "tools.h"

using namespace Rcpp;
//using namespace arma;
using namespace std;
//using namespace Eigen;
//using namespace Numer;

// USE A WAY TO QUICLY COMPUTE THE DET AND THE INVERSE 
// A Structural Model for the Coevolution of Networks and Behavior

// without contextual effects

// [[Rcpp::export]]
List peerMCMCnoc_none(const List& y,
                  const List& V,   // is Xone
                  List& Gnorm,
                  List& prior,
                  List ListIndex,
                  const int& M,
                  const IntegerVector& N,
                  const int& kbeta,
                  const arma::vec& theta0,
                  const arma::mat& invsigmatheta,
                  const double& zeta0,
                  const double& invsigma2zeta,
                  const double& a, 
                  const double& b,
                  const arma::vec& parms0,
                  const int& iteration,
                  const double& target,
                  const double& jumpmin,
                  const double& jumpmax,
                  const double& c,
                  const int& progress){
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  int sumN = sum(N);
  
  int kv = kbeta; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  
  // Other parameters
  double jumpzeta = 1;
  double zetaaccept = 0;
  List Vtheta(M), A(M), Ay(M);
  double sumlogdetA = 0.0;
  for(int m(0); m<M; ++m){
    //Gnorm
    int Nm = N(m);
    arma::mat Gm = Gnorm[m];
    
    // V
    arma::mat Vm = V[m];
    arma::vec Vmtheta = Vm*theta;
    // Vtheta
    Vtheta[m] = Vmtheta;
    
    // A
    arma::mat Am = arma::eye(Nm,Nm) - alpha*Gm;
    A[m] = Am;
    
    // Ay
    arma::vec ym = y[m];
    Ay[m] = Am*ym;
    
    // sumlogdetA
    double valm, signm;
    arma::log_det(valm, signm, Am);
    sumlogdetA += valm + log(signm);
  }
  
  
  
  //Save 
  arma::mat saveparms(kv+2,iteration);
  NumericVector parmscpp;
  // loop
  
  if (progress == 0 ){
    
    for(int t(0); t<iteration; ++t){
      // Update G
      updGnormnoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
    }
  }
  if (progress == 2 ){
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      
      
      // Update G
      updGnormnoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      parmscpp             = wrap(parms);
      parmscpp.attr("dim") = R_NilValue;
      Rcpp::print(parmscpp);
      //Rcpp::Rcout<<"************************"<<std::endl;
      Rprintf("************************ \n");
    }}
  if (progress == 1 ){
    Progress prgcpp(iteration, true);
    for(int t(0); t<iteration; ++t){
      prgcpp.increment();

      // Update G
      updGnormnoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
    }
  }
  
  
  
  NumericMatrix output = wrap(saveparms.t());
  
  return List::create(Named("posterior")  = output, 
                      Named("acceptance") = zetaaccept/iteration,
                      Named("G")          = Gnorm);
}


// [[Rcpp::export]]
List peerMCMCblocknoc_none(const List& y,
                       const List& V,   // is Xone
                       List& Gnorm,
                       List& prior,
                       List ListIndex,
                       const int& M,
                       const IntegerVector& N,
                       const int& kbeta,
                       const arma::vec& theta0,
                       const arma::mat& invsigmatheta,
                       const double& zeta0,
                       const double& invsigma2zeta,
                       const double& a, 
                       const double& b,
                       const arma::vec& parms0,
                       const int& iteration,
                       const double& target,
                       const double& jumpmin,
                       const double& jumpmax,
                       const double& c,
                       const int& progress,
                       const int& nupmax){
  
  
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  int sumN = sum(N);
  
  int kv = kbeta; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  
  // Other parameters
  double jumpzeta = 1;
  double zetaaccept = 0;
  List Vtheta(M), A(M), Ay(M);
  double sumlogdetA = 0.0;
  for(int m(0); m<M; ++m){
    //Gnorm
    int Nm = N(m);
    arma::mat Gm = Gnorm[m];
    
    // V
    arma::mat Vm = V[m];
    arma::vec Vmtheta = Vm*theta;
    // Vtheta
    Vtheta[m] = Vmtheta;
    
    // A
    arma::mat Am = arma::eye(Nm,Nm) - alpha*Gm;
    A[m] = Am;
    
    // Ay
    arma::vec ym = y[m];
    Ay[m] = Am*ym;
    
    // sumlogdetA
    double valm, signm;
    arma::log_det(valm, signm, Am);
    sumlogdetA += valm + log(signm);
  }
  
  
  //Save 
  arma::mat saveparms(kv+2,iteration);
  NumericVector parmscpp;
  // loop
  if (progress == 0 ){
    
    for(int t(0); t<iteration; ++t){
      
      // Update G
      updGnormblocknoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2, nupmax);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
    }
  }
  if (progress == 2){
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);

      // Update G
      updGnormblocknoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2, nupmax);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      parmscpp             = wrap(parms);
      parmscpp.attr("dim") = R_NilValue;
      Rcpp::print(parmscpp);
      //Rcpp::Rcout<<"************************"<<std::endl;
      Rprintf("************************ \n");
    }
  }
  if (progress == 1){
    Progress prgcpp(iteration, true);
    for(int t(0); t<iteration; ++t){
      prgcpp.increment();

      // Update G
      updGnormblocknoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2, nupmax);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
    }
  }
  
  NumericMatrix output = wrap(saveparms.t());
  
  return List::create(Named("posterior")  = output, 
                      Named("acceptance") = zetaaccept/iteration,
                      Named("G")          = Gnorm);
}

//////////////////////////////// with contextual effects


// [[Rcpp::export]]
List peerMCMC_none(const List& y,
               const List& X,
               const List& Xone,
               List& Gnorm,
               List& prior,
               List ListIndex,
               const int& M,
               const IntegerVector& N,
               const int& kbeta,
               const int& kgamma,
               const arma::vec& theta0,
               const arma::mat& invsigmatheta,
               const double& zeta0,
               const double& invsigma2zeta,
               const double& a, 
               const double& b,
               const arma::vec& parms0,
               const int& iteration,
               const double& target,
               const double& jumpmin,
               const double& jumpmax,
               const double& c,
               const int& progress){
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  int sumN = sum(N);
  
  int kv = kbeta + kgamma; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  
  // Other parameters
  double jumpzeta = 1;
  double zetaaccept = 0;
  List V(M), Vtheta(M), Xb(M), Xgamma(M), A(M), Ay(M);
  double sumlogdetA = 0.0;
  for(int m(0); m<M; ++m){
    //Gnorm
    int Nm = N(m);
    arma::mat Gm = Gnorm[m];
    
    // V
    arma::mat Xonem = Xone[m];
    arma::mat Xm = X[m];
    arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
    arma::vec Vmtheta = Vm*theta;
    arma::vec Xmb = Xonem*theta.head(kbeta);
    arma::vec Xmgamma = Xm*theta.tail(kgamma);
    V[m] = Vm;
    
    // Vtheta
    Vtheta[m] = Vmtheta;
    
    //Xb
    Xb[m] = Xmb;
    
    //Xgamma
    Xgamma[m] = Xmgamma;
    
    // A
    arma::mat Am = arma::eye(Nm,Nm) - alpha*Gm;
    A[m] = Am;
    
    // Ay
    arma::vec ym = y[m];
    Ay[m] = Am*ym;
    
    // sumlogdetA
    double valm, signm;
    arma::log_det(valm, signm, Am);
    sumlogdetA += valm + log(signm);
  }
  
  
  
  //Save 
  arma::mat saveparms(kv+2,iteration);
  NumericVector parmscpp;
  // loop
  
  if (progress == 0 ){
    
    for(int t(0); t<iteration; ++t){
      // Update G
      updGnorm (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
        Vtheta[m] = Vmtheta;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
    }
  }
  if (progress == 2 ){
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      
      // Update G
      updGnorm (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
        Vtheta[m] = Vmtheta;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      parmscpp             = wrap(parms);
      parmscpp.attr("dim") = R_NilValue;
      Rcpp::print(parmscpp);
      //Rcpp::Rcout<<"************************"<<std::endl;
      Rprintf("************************ \n");
    }}
  if (progress == 1 ){
    Progress prgcpp(iteration, true);
    for(int t(0); t<iteration; ++t){
      prgcpp.increment();
      // Update G
      updGnorm (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
        Vtheta[m] = Vmtheta;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
    }
  }
  
  
  NumericMatrix output = wrap(saveparms.t());
  
  return List::create(Named("posterior")  = output, 
                      Named("acceptance") = zetaaccept/iteration,
                      Named("G")          = Gnorm);
}


// [[Rcpp::export]]
List peerMCMCblock_none(const List& y,
                    const List& X,
                    const List& Xone,
                    List& Gnorm,
                    List& prior,
                    List ListIndex,
                    const int& M,
                    const IntegerVector& N,
                    const int& kbeta,
                    const int& kgamma,
                    const arma::vec& theta0,
                    const arma::mat& invsigmatheta,
                    const double& zeta0,
                    const double& invsigma2zeta,
                    const double& a, 
                    const double& b,
                    const arma::vec& parms0,
                    const int& iteration,
                    const double& target,
                    const double& jumpmin,
                    const double& jumpmax,
                    const double& c, 
                    const int& progress,
                    const int& nupmax){
  
  
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  int sumN = sum(N);
  
  int kv = kbeta + kgamma; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  
  // Other parameters
  double jumpzeta = 1;
  double zetaaccept = 0;
  List V(M), Vtheta(M), Xb(M), Xgamma(M), A(M), Ay(M);
  double sumlogdetA = 0.0;
  for(int m(0); m<M; ++m){
    //Gnorm
    int Nm = N(m);
    arma::mat Gm = Gnorm[m];
    
    // V
    arma::mat Xonem = Xone[m];
    arma::mat Xm = X[m];
    arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
    arma::vec Vmtheta = Vm*theta;
    arma::vec Xmb = Xonem*theta.head(kbeta);
    arma::vec Xmgamma = Xm*theta.tail(kgamma);
    V[m] = Vm;
    
    // Vtheta
    Vtheta[m] = Vmtheta;
    
    //Xb
    Xb[m] = Xmb;
    
    //Xgamma
    Xgamma[m] = Xmgamma;
    
    // A
    arma::mat Am = arma::eye(Nm,Nm) - alpha*Gm;
    A[m] = Am;
    
    // Ay
    arma::vec ym = y[m];
    Ay[m] = Am*ym;
    
    // sumlogdetA
    double valm, signm;
    arma::log_det(valm, signm, Am);
    sumlogdetA += valm + log(signm);
  }
  
  
  //Save 
  arma::mat saveparms(kv+2,iteration);
  NumericVector parmscpp;
  // loop
  if (progress == 0 ){
    
    for(int t(0); t<iteration; ++t){
      // Update G
      updGnormblock (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2, nupmax);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
        Vtheta[m] = Vmtheta;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
    }
  }
  if (progress == 2 ) {
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      // Update G
      updGnormblock (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2, nupmax);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
        Vtheta[m] = Vmtheta;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      parmscpp             = wrap(parms);
      parmscpp.attr("dim") = R_NilValue;
      Rcpp::print(parmscpp);
      //Rcpp::Rcout<<"************************"<<std::endl;
      Rprintf("************************ \n");
    }
  }
  if (progress == 1 ){
    Progress prgcpp(iteration, true);
    for(int t(0); t<iteration; ++t){
      prgcpp.increment();
      // Update G
      updGnormblock (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2, nupmax);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
        Vtheta[m] = Vmtheta;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target)/pow(t+1,c);
      fsetjump_d(jumpzeta, jumpmin, jumpmax);
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
    }
  }
  
  
  NumericMatrix output = wrap(saveparms.t());
  
  return List::create(Named("posterior")  = output, 
                      Named("acceptance") = zetaaccept/iteration,
                      Named("G")          = Gnorm);
}