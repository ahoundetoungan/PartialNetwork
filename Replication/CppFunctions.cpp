// [[Rcpp::depends(RcppEigen, RcppArmadillo, RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>

using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;


// This function computes the likelihood for a consored Poisson regression
//[[Rcpp::export]]
double f_rcpoisson(const arma::vec& beta, 
                   const arma::vec& y,
                   const arma::mat& X,
                   const Rcpp::LogicalVector& censure,
                   const arma::vec& lcensure,
                   const arma::vec& rcensure){
  // cout<<"beta: "<<beta.t()<<endl;
  arma::vec eXb = exp(X*beta);
  int n         = y.n_elem;
  double llh = 0;
  for(int i(0); i < n; ++ i){
    if(censure(i)){
      double l2 = R::ppois(rcensure(i), eXb(i), true, true);
      double l1 = R::ppois(lcensure(i), eXb(i), true, true);
      llh      += (l1 + log(exp(l2 - l1) - 1));
    } else {
      llh      += R::dpois(y(i), eXb(i), true);
    }
  }
  return -llh;
}

// This function computes the expectation of y when y follows a consored Poisson distribution
//[[Rcpp::export]]
Rcpp::NumericVector expy(const arma::vec& beta, 
                         const arma::mat& X,
                         const arma::vec& lcensure,
                         const arma::vec& rcensure){
  // cout<<"beta: "<<beta.t()<<endl;
  arma::vec Xb    = X*beta;
  arma::vec eXb   = exp(Xb);
  int n           = X.n_rows;
  Rcpp::NumericVector out(n);
  for(int i(0); i < n; ++ i){
    double lpa    = R::ppois(lcensure(i), eXb(i), true, true);
    double lpamo  = R::ppois(lcensure(i) - 1, eXb(i), true, true);
    double lpb    = R::ppois(rcensure(i), eXb(i), true, true);
    double lpbmo  = R::ppois(rcensure(i) - 1, eXb(i), true, true);
    out(i)        = Xb(i) + lpamo - lpa + log(exp(lpbmo - lpamo) - 1) - log(exp(lpb - lpa) - 1);
  }
  out.attr("dim") = R_NilValue;
  return exp(out);
}




// This is a class to estimate a logit model with misclassification
// Type 1: Only sensitivity is model, specificity = 1
class LogitMisclassified1: public MFuncGrad
{
private:
  const Eigen::MatrixXd& X;
  const Eigen::ArrayXd& y;
  const int& Kx;
public:
  LogitMisclassified1(const Eigen::MatrixXd& X_, 
                     const Eigen::ArrayXd& y_, 
                     const int& Kx_) 
    : X(X_), y(y_), Kx(Kx_) {}
  Eigen::VectorXd Grad;
  double f_grad(Constvec& theta, Refvec grad) {
    double p11(1.0/(1 + exp(-theta(0)))); // sensitivity
    Eigen::VectorXd beta(theta.segment(1, Kx));
    
    Eigen::ArrayXd xb(X * beta);
    Eigen::ArrayXd exb(xb.exp());
    
    Eigen::ArrayXd tp1(p11*exb);
    Eigen::ArrayXd tp2(1 + (1 - p11)*exb);
    Eigen::ArrayXd tp3(1 + exb);
    
    // likelihood
    double llh = (y*tp1.log() + (1 - y)*tp2.log() - tp3.log()).sum();
    
    // Gradient
    Grad    = Eigen::VectorXd::Zero(Kx + 1);
    Grad(0) = ((y/tp1 - (1 - y)/tp2)*exb).sum()*exp(theta(1))/pow(1 + exp(theta(1)), 2);
    Grad.segment(1, Kx) = X.transpose()*((p11*y/tp1 + (1 - p11)*(1 - y)/tp2 - 1/tp3)*exb).matrix();
    grad = -Grad;
    return -llh;
  }
};

// Type 2: Only specificity is model, sensitivity = 1
class LogitMisclassified2: public MFuncGrad
{
private:
  const Eigen::MatrixXd& X;
  const Eigen::ArrayXd& y;
  const int& Kx;
public:
  LogitMisclassified2(const Eigen::MatrixXd& X_, 
                     const Eigen::ArrayXd& y_, 
                     const int& Kx_) 
    : X(X_), y(y_), Kx(Kx_) {}
  Eigen::VectorXd Grad;
  double f_grad(Constvec& theta, Refvec grad) {
    double p00(1.0/(1 + exp(-theta(0)))); // specificity
    Eigen::VectorXd beta(theta.segment(1, Kx));
    
    Eigen::ArrayXd xb(X * beta);
    Eigen::ArrayXd exb(xb.exp());
    
    Eigen::ArrayXd tp1(1 - p00 + exb);
    Eigen::ArrayXd tp2(p00*Eigen::ArrayXd::Ones(exb.size()));
    Eigen::ArrayXd tp3(1 + exb);
    
    // likelihood
    double llh = (y*tp1.log() + (1 - y)*tp2.log() - tp3.log()).sum();
    
    // Gradient
    Grad    = Eigen::VectorXd::Zero(Kx + 1);
    Grad(0) = ((1 - y)/tp2 - y/tp1).sum()*exp(theta(0))/pow(1 + exp(theta(0)), 2);
    Grad.segment(1, Kx) = X.transpose()*((y/tp1 - 1/tp3)*exb).matrix();
    grad = -Grad;
    return -llh;
  }
};

// Type 3: Both specificity and sensitivity are model
class LogitMisclassified3: public MFuncGrad
{
private:
  const Eigen::MatrixXd& X;
  const Eigen::ArrayXd& y;
  const int& Kx;
public:
  LogitMisclassified3(const Eigen::MatrixXd& X_, 
                     const Eigen::ArrayXd& y_, 
                     const int& Kx_) 
    : X(X_), y(y_), Kx(Kx_) {}
  Eigen::VectorXd Grad;
  double f_grad(Constvec& theta, Refvec grad) {
    double p00(1.0/(1 + exp(-theta(0)))); // specificity
    double p11(1.0/(1 + exp(-theta(1)))); // sensitivity
    Eigen::VectorXd beta(theta.segment(2, Kx));
    
    Eigen::ArrayXd xb(X * beta);
    Eigen::ArrayXd exb(xb.exp());
    
    Eigen::ArrayXd tp1(1 - p00 + p11*exb);
    Eigen::ArrayXd tp2(p00 + (1 - p11)*exb);
    Eigen::ArrayXd tp3(1 + exb);
    
    // likelihood
    double llh = (y*tp1.log() + (1 - y)*tp2.log() - tp3.log()).sum();
    
    // Gradient
    Grad    = Eigen::VectorXd::Zero(Kx + 2);
    Grad(0) = ((1 - y)/tp2 - y/tp1).sum()*exp(theta(0))/pow(1 + exp(theta(0)), 2);
    Grad(1) = ((y/tp1 - (1 - y)/tp2)*exb).sum()*exp(theta(1))/pow(1 + exp(theta(1)), 2);
    Grad.segment(2, Kx) = X.transpose()*((p11*y/tp1 + (1 - p11)*(1 - y)/tp2 - 1/tp3)*exb).matrix();
    grad = -Grad;
    return -llh;
  }
};

// This function estimates a logit model with misclassification
// [[Rcpp::export]]
Rcpp::List flogistmisclassify(const Eigen::VectorXd theta,
                              const Eigen::ArrayXd& y,  
                              const Eigen::MatrixXd& X, 
                              const int type  = 3, // 1 is only sensitivity, 2 only specificity, 3 is both
                              const int maxit = 1e9, 
                              const double& eps_f = 1e-12, 
                              const double& eps_g = 1e-12){
  int Kx(X.cols()), status;
  double fopt;
  Eigen::VectorXd theta_(theta), Grad;
  
  // if (type == 1) {
  //   if(theta_.size() != (Kx + 1)) Rcpp::stop("The size of `theta` does not match `type`.");
  //   LogitMisclassified1 nll(X, y, Kx);
  //   theta_(0) = log(theta_(0)/(1 - theta_(0)));
  //   status    = optim_lbfgs(nll, theta_, fopt, maxit, eps_f, eps_g);
  //   theta_(0) = 1/(1 + exp(-theta_(0)));
  //   Grad      = nll.Grad;
  // } else if (type == 2) {
  //   if(theta_.size() != (Kx + 1)) Rcpp::stop("The size of `theta` does not match `type`.");
  //   LogitMisclassified2 nll(X, y, Kx);
  //   theta_(0) = log(theta_(0)/(1 - theta_(0)));
  //   status    = optim_lbfgs(nll, theta_, fopt, maxit, eps_f, eps_g);
  //   theta_(0) = 1/(1 + exp(-theta_(0)));
  //   Grad      = nll.Grad;
  // } else if (type == 3) {
  //   if(theta_.size() != (Kx + 1)) Rcpp::stop("The size of `theta` does not match `type`.");
  //   LogitMisclassified3 nll(X, y, Kx);
  //   theta_(0) = log(theta_(0)/(1 - theta_(0)));
  //   theta_(1) = log(theta_(1)/(1 - theta_(1)));
  //   status    = optim_lbfgs(nll, theta_, fopt, maxit, eps_f, eps_g);
  //   theta_(0) = 1/(1 + exp(-theta_(0)));
  //   theta_(1) = 1/(1 + exp(-theta_(1)));
  //   Grad      = nll.Grad;
  // }
  
  return Rcpp::List::create(
    Rcpp::Named("estimate") = theta_,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = Grad,
    Rcpp::Named("status")   = status);
}
