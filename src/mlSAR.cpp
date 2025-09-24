// [[Rcpp::depends(RcppArmadillo, RcppNumerical)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppNumerical.h>

using namespace Rcpp;
//using namespace arma;
using namespace std;
//using namespace Eigen;
//using namespace Numer;

class Peer: public Numer::MFuncGrad
{
private:
  List& listG;
  const Rcpp::IntegerVector& N;
  const int& M;
  const List& y;
  const int& kz; //  number of parameters excepted sigma2 and alpha
public:
  Peer(List& listG_, const Rcpp::IntegerVector& N_, const int& M_, const List& y_, const int& kz_) : 
  listG(listG_), N(N_), M(M_), y(y_), kz(kz_){}
  
  double alpha;
  arma::vec beta;
  double sigma2;
  arma::vec dbeta;
  double dsigma2;
  double Grad;
  
  
  double f_grad(Numer::Constvec& parms, Numer::Refvec grad)
  {
    Eigen::VectorXd parms0 = parms;  //make a copy
    arma::vec alphaAR = arma::vec(parms0 .data(), parms0 .size(), false, false); //converte into arma vec
    
    alpha = exp(alphaAR(0))/(1+exp(alphaAR(0)));
    
    List Gy = listG(0), Z = listG(1), EI = listG(2);
    
    arma::mat sZZ(kz,kz,arma::fill::zeros);
    arma::vec sZAY(kz,arma::fill::zeros);
    arma::vec sZGY(kz,arma::fill::zeros);
    
    double sumN = sum(N);
    //compute beta as [beta,gamma]
    
    for(int m(0); m<M; ++m){
      arma::vec Gym = Gy[m], ym=y[m];
      
      arma::mat Zm = Z[m];
      
      sZZ += Zm.t()*Zm;
      sZAY += Zm.t()*(ym-alpha*Gym);
      sZGY += Zm.t()*(-Gym);
    }
    
    beta = arma::solve(sZZ, sZAY);
    dbeta = arma::solve(sZZ, sZGY);
    
    //compute sigma2
    double ee=0.0;
    
    arma::vec sA(1,arma::fill::zeros);
    arma::rowvec sB(kz,arma::fill::zeros), sD(kz,arma::fill::zeros);
    arma::mat sC(kz,kz,arma::fill::zeros);
    
    for(int m(0); m<M; ++m){
      arma::vec Gym = Gy[m], ym=y[m];
      
      arma::mat Zm = Z[m];
      
      sA += (arma::trans(alpha*Gym - ym)*Gym);
      sB += Gym.t()*Zm;
      sC += Zm.t()*Zm;
      sD += ym.t()*Zm;
      
      arma::vec em = ym - alpha*Gym - Zm*beta;
      
      ee += dot(em,em);
    }
    
    sigma2 = ee/sumN;
    dsigma2 = arma::as_scalar(sA + sB*(beta + alpha*dbeta) + (beta.t()*sC - sD)*dbeta);
    
    double llh=0.0;
    double E=0.0;
    
    for(int m(0); m<M; ++m){
      arma::cx_vec EIGENm = EI[m];
      arma::cx_vec tempeig = alpha*EIGENm - 1;
      
      llh += log(abs(real(arma::prod(tempeig))));
      E +=  real(sum(EIGENm/tempeig));
    }
    
    llh = 0.5*sumN*(log(sigma2)+1) - llh;
    
    Grad = (dsigma2/sigma2 - E)*alpha/((1+exp(alphaAR(0))));
    
    grad[0]= Grad/sumN;
    
    
    return llh/sumN;
  }
};


List estimatepeerML(List& listG,
                    const Rcpp::IntegerVector& N,
                    const int& M,
                    const List& y,
                    const int& kz) {
  // Negative log likelihood
  Peer f(listG, N, M, y, kz);
  
  // Initial guess
  Eigen::VectorXd transalpha(1);
  transalpha(0)=0;
  
  double fopt;
  int status = optim_lbfgs(f, transalpha, fopt);
  
  //End Temporary
  
  return Rcpp::List::create(
    Rcpp::Named("alpha") = f.alpha,
    Rcpp::Named("beta") = f.beta,
    Rcpp::Named("sigma2") = f.sigma2,
    Rcpp::Named("status") = status);
}


// Starting point
// [[Rcpp::export]] 
arma::vec sartpoint (List& Gnorm,
                     const int& M,
                     const IntegerVector& N,
                     const int& kbeta,
                     const int& kgamma,
                     const List& y,
                     const List& X, 
                     const List& Xone) {
  
  
  int kz = kbeta + kgamma; // number of exogenous variables
  
  
  List listG(3);
  
  //Save results
  
  arma::vec theta(kz+2);
  List Gy(M), Z(M), EIOGGX(M);
  
  for(int m(0); m<M; ++m){
    arma::mat Gm = Gnorm(m);
    arma::mat Xm = Xone(m);
    arma::mat Xmshed = X(m);
    arma::vec ym = y(m);
    
    // GY
    Gy(m) = Gm*ym;
    
    // GX
    Z(m) = arma::join_rows(Xm, Gm*Xmshed);
    
    // Eigen
    arma::cx_vec EIGENm;
    arma::eig_gen(EIGENm, Gm);
    EIOGGX(m) = EIGENm;
  }
  
  listG(0) = Gy;
  listG(1) = Z;
  listG(2) = EIOGGX;
  
  List MLestim=estimatepeerML(listG, N, M, y, kz);
  
  double alphaML = MLestim(0);
  arma::vec thetaML = MLestim(1);
  double seML = MLestim(2);
  
  theta.head(kz) = thetaML;
  theta(kz)=alphaML;
  theta(kz+1)=seML;
  
  return theta;
}


// [[Rcpp::export]] 
arma::vec sartpointnoc (List& Gnorm,
                        const int& M,
                        const IntegerVector& N,
                        const int& kbeta,
                        const List& y,
                        const List& Xone) {
  
  
  int kz = kbeta; // number of exogenous variables
  
  
  List listG(3);
  
  //Save results
  
  arma::vec theta(kz+2);
  List Gy(M), EIOGGX(M);
  
  for(int m(0); m<M; ++m){
    arma::mat Gm = Gnorm(m);
    arma::vec ym = y(m);
    
    // GY
    Gy(m) = Gm*ym;
    
    // Eigen
    arma::cx_vec EIGENm;
    arma::eig_gen(EIGENm, Gm);
    EIOGGX(m) = EIGENm;
  }
  
  listG(0) = Gy;
  listG(1) = Xone;
  listG(2) = EIOGGX;
  
  List MLestim=estimatepeerML(listG, N, M, y, kz);
  
  double alphaML = MLestim(0);
  arma::vec thetaML = MLestim(1);
  double seML = MLestim(2);
  
  theta.head(kz) = thetaML;
  theta(kz)=alphaML;
  theta(kz+1)=seML;
  
  return theta;
}
