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
List peerMCMCnoc_pl(const List& y,
                    const List& V,   // is Xone
                    List& Gnorm,
                    List& G0obs,
                    List& prior,
                    List& ListIndex,
                    const int& M,
                    const IntegerVector& N,
                    const int& kbeta,
                    const arma::vec& theta0,
                    const arma::mat& invsigmatheta,
                    const double& zeta0,
                    const double& invsigma2zeta,
                    const double& a, 
                    const double& b,
                    const Rcpp::NumericVector weight,
                    const arma::mat& dZ,
                    const arma::vec& murho,
                    const arma::mat& Vrho,
                    const int& Krho,
                    Eigen::VectorXd& lFdZrho1,
                    Eigen::VectorXd& lFdZrho0,
                    const arma::vec& parms0,
                    const int& iteration,
                    const arma::vec& target,
                    const arma::vec& jumpmin,
                    const arma::vec& jumpmax,
                    const double& c,
                    const int& progress,
                    const int& type,
                    const bool& Afixed,
                    const Eigen::ArrayXd& G0obsvec){
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  int sumN = sum(N);
  
  int kv = kbeta; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  arma::vec rho = murho;
  arma::mat jumprho = Vrho;
  arma::mat iVrho   = arma::inv(Vrho);
  
  // initialize network model's parameters
  Eigen::VectorXd lFdZrhoE1  = lFdZrho1;
  Eigen::VectorXd lFdZrhoE0  = lFdZrho0;
  
  
  // Other parameters
  double jumpzeta   = 1;
  double zetaaccept = 0;
  double djumprho   = 1;
  double rhoaccept  = 0;
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
  arma::mat saverho(Krho,iteration);
  NumericVector parmscpp;
  // loop
  
  if (progress == 0 ){
    
    for(int t(0); t<iteration; ++t){
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
      // Update G
      updGnormnoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      
      saverho.col(t)   = rho;
    }
  }
  if (progress == 2 ){
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
      // Update G
      updGnormnoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
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
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
      // Update G
      updGnormnoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
    }
  }
  
  
  
  List output              = List::create(Named("theta") = wrap(saveparms.t()),
                                          Named("rho")   = wrap(saverho.t()));
  NumericVector acceptance = NumericVector::create(_["alpha"] = zetaaccept/iteration,
                                                   _["rho"]          = rhoaccept/iteration);
  
  return List::create(Named("posterior")  = output, 
                      Named("acceptance") = acceptance,
                      Named("G")          = Gnorm);
}


// [[Rcpp::export]]
List peerMCMCblocknoc_pl(const List& y,
                         const List& V,   // is Xone
                         List& Gnorm,
                         List& G0obs,
                         List& prior,
                         List& ListIndex,
                         const int& M,
                         const IntegerVector& N,
                         const int& kbeta,
                         const arma::vec& theta0,
                         const arma::mat& invsigmatheta,
                         const double& zeta0,
                         const double& invsigma2zeta,
                         const double& a, 
                         const double& b,
                         const Rcpp::NumericVector weight,
                         const arma::mat& dZ,
                         const arma::vec& murho,
                         const arma::mat& Vrho,
                         const int& Krho,
                         Eigen::VectorXd& lFdZrho1,
                         Eigen::VectorXd& lFdZrho0,
                         const arma::vec& parms0,
                         const int& iteration,
                         const arma::vec& target,
                         const arma::vec& jumpmin,
                         const arma::vec& jumpmax,
                         const double& c,
                         const int& progress,
                         const int& nupmax,
                         const int& type,
                         const bool& Afixed,
                         const Eigen::ArrayXd& G0obsvec){
  
  
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  int sumN = sum(N);
  
  int kv = kbeta; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  arma::vec rho = murho;
  arma::mat jumprho = Vrho;
  arma::mat iVrho   = arma::inv(Vrho);
  
  // initialize network model's parameters
  Eigen::VectorXd lFdZrhoE1  = lFdZrho1;
  Eigen::VectorXd lFdZrhoE0  = lFdZrho0;
  
  
  // Other parameters
  double jumpzeta   = 1;
  double zetaaccept = 0;
  double djumprho   = 1;
  double rhoaccept  = 0;
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
  arma::mat saverho(Krho,iteration);
  NumericVector parmscpp;
  // loop
  if (progress == 0 ){
    
    for(int t(0); t<iteration; ++t){
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
      // Update G
      updGnormblocknoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2, nupmax);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
    }
  }
  if (progress == 2){
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
      // Update G
      updGnormblocknoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2, nupmax);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      // save output
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
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
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
      // Update G
      updGnormblocknoc (Gnorm, prior,ListIndex, N, M, y, A, Ay, Vtheta, alpha, sigma2, nupmax);
      
      // Update theta
      updthetanoc (theta, Vtheta, sigma2, kv, Ay, V, invsigmathetatheta0,
                   invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
    }
  }
  
  List output              = List::create(Named("theta") = wrap(saveparms.t()),
                                          Named("rho")   = wrap(saverho.t()));
  NumericVector acceptance = NumericVector::create(_["alpha"] = zetaaccept/iteration,
                                                   _["rho"]          = rhoaccept/iteration);
  
  return List::create(Named("posterior")  = output, 
                      Named("acceptance") = acceptance,
                      Named("G")          = Gnorm);
}
//////////////////////////////// with contextual effects


// [[Rcpp::export]]
List peerMCMC_pl(const List& y,
                 const List& X,
                 const List& Xone,
                 List& Gnorm,
                 List& G0obs,
                 List& prior,
                 List& ListIndex,
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
                 const arma::mat& dZ,
                 const Rcpp::NumericVector weight,
                 const arma::vec& murho,
                 const arma::mat& Vrho,
                 const int& Krho,
                 Eigen::VectorXd& lFdZrho1,
                 Eigen::VectorXd& lFdZrho0,
                 const arma::vec& parms0,
                 const int& iteration,
                 const arma::vec& target,
                 const arma::vec& jumpmin,
                 const arma::vec& jumpmax,
                 const double& c,
                 const int& progress,
                 const int& type,
                 const bool& Afixed,
                 const Eigen::ArrayXd& G0obsvec){
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  int sumN = sum(N);
  
  int kv = kbeta + kgamma; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  arma::vec rho = murho;
  arma::mat jumprho = Vrho;
  arma::mat iVrho   = arma::inv(Vrho);
  
  // initialize network model's parameters
  Eigen::VectorXd lFdZrhoE1  = lFdZrho1;
  Eigen::VectorXd lFdZrhoE0  = lFdZrho0;
  
  
  // Other parameters
  double jumpzeta   = 1;
  double zetaaccept = 0;
  double djumprho   = 1;
  double rhoaccept  = 0;
  
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
  arma::mat saverho(Krho,iteration);
  NumericVector parmscpp;
  // loop
  
  if (progress == 0 ){
    
    for(int t(0); t<iteration; ++t){
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
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
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
    }
  }
  if (progress == 2 ){
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
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
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
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
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
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
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
    }
  }
  
  
  List output              = List::create(Named("theta") = wrap(saveparms.t()),
                                          Named("rho")   = wrap(saverho.t()));
  NumericVector acceptance = NumericVector::create(_["alpha"] = zetaaccept/iteration,
                                                   _["rho"]          = rhoaccept/iteration);
  
  return List::create(Named("posterior")  = output, 
                      Named("acceptance") = acceptance,
                      Named("G")          = Gnorm);
}


// [[Rcpp::export]]
List peerMCMCblock_pl(const List& y,
                      const List& X,
                      const List& Xone,
                      List& Gnorm,
                      List& G0obs,
                      List& prior,
                      List& ListIndex,
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
                      const Rcpp::NumericVector weight,
                      const arma::mat& dZ,
                      const arma::vec murho,
                      const arma::mat& Vrho,
                      const int& Krho,
                      Eigen::VectorXd& lFdZrho1,
                      Eigen::VectorXd& lFdZrho0,
                      const arma::vec& parms0,
                      const int& iteration,
                      const arma::vec& target,
                      const arma::vec& jumpmin,
                      const arma::vec& jumpmax,
                      const double& c, 
                      const int& progress,
                      const int& nupmax,
                      const int& type,
                      const bool& Afixed,
                      const Eigen::ArrayXd& G0obsvec){
  
  
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  int sumN = sum(N);
  
  int kv = kbeta + kgamma; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  arma::vec rho = murho;
  arma::mat jumprho = Vrho;
  arma::mat iVrho   = arma::inv(Vrho);
  
  // initialize network model's parameters
  Eigen::VectorXd lFdZrhoE1  = lFdZrho1;
  Eigen::VectorXd lFdZrhoE0  = lFdZrho0;
  
  
  // Other parameters
  double jumpzeta   = 1;
  double zetaaccept = 0;
  double djumprho   = 1;
  double rhoaccept  = 0;
  
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
  arma::mat saverho(Krho,iteration);
  NumericVector parmscpp;
  // loop
  if (progress == 0 ){
    
    for(int t(0); t<iteration; ++t){
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
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
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
    }
  }
  if (progress == 2 ) {
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
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
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
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
      //Update prior
      updrhopl(Gnorm, prior, G0obs, ListIndex, rho, lFdZrhoE1, lFdZrhoE0, weight, dZ, murho, iVrho, jumprho,
               Krho,  N,  M,  rhoaccept, type, Afixed, G0obsvec);
      djumprho += (rhoaccept/(t + 1) - target(1))/pow(t + 1,c);
      fsetjump_dm(djumprho, jumpmin(1), jumpmax(1), jumprho, Vrho);
      
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
      jumpzeta += (zetaaccept/(t + 1) - target(0))/pow(t + 1,c);
      fsetjump_d(jumpzeta, jumpmin(0), jumpmax(0));
      
      // save output
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      saverho.col(t)   = rho;
    }
  }
  
  
  List output              = List::create(Named("theta") = wrap(saveparms.t()),
                                          Named("rho")   = wrap(saverho.t()));
  NumericVector acceptance = NumericVector::create(_["alpha"] = zetaaccept/iteration,
                                                   _["rho"]          = rhoaccept/iteration);
  
  return List::create(Named("posterior")  = output, 
                      Named("acceptance") = acceptance,
                      Named("G")          = Gnorm);
}