#ifndef PACKAGECPdetect_TOO_H
#define PACKAGECPdetect_TOO_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>

arma::vec Fmvnorm(const double& dim, arma::vec u, arma::mat sigma);

Eigen::MatrixXd invmodij (const Eigen::MatrixXd& invM,
                          const unsigned int& i,
                          const unsigned int& j,
                          const double& eps_);
Eigen::MatrixXd invmodijk (const Eigen::MatrixXd& invM,
                           const unsigned int& i,
                           const unsigned int& j,
                           const unsigned int& k,
                           const double& eps_1,
                           const double& eps_2);
double detmodij (const double& detM,
                 const Eigen::MatrixXd& invM,
                 const unsigned int& i,
                 const unsigned int& j,
                 const double& eps_) ;
double detmodijk (const double& detM,
                  const Eigen::MatrixXd& invM,
                  const unsigned int& i,
                  const unsigned int& j,
                  const unsigned int& k,
                  const double& eps_1,
                  const double& eps_2);
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
arma::mat possentries(int& nupdate, const int& pow_nupdate);
void updselel (
        Eigen::MatrixXd& Gmnorm,
        double& prior_blockl,
        const Eigen::MatrixXd priorm,
        const int& i,
        const arma::vec& index_col,
        const arma::rowvec& newval,
        const int& nupdate
);
arma::vec cBlock(
        const int& nupmax, 
        const int& NJ, 
        int& NJeff);

Rcpp::List frVtoM(const Eigen::VectorXd& u,
                  const Rcpp::IntegerVector& N,
                  const double& M);

Rcpp::List frVtoMarma(const arma::vec& u,
                const Rcpp::IntegerVector& N,
                const double& M);

Eigen::VectorXd frMtoV(Rcpp::List& u,
                       const Rcpp::IntegerVector& N,
                       const double& M);

Rcpp::List fListIndex(Rcpp::List& prior,
                      Rcpp::List& G0obs,
                      const int& M,
                      const Rcpp::IntegerVector& N);
void fsetjump_dm(double& jump, const double& jumpmin,
                const double& jumpmax, arma::mat& jumpl, const arma::mat& Vparms);
void fsetjump_vl(arma::vec& jump, const double& jumpmin,
            const double& jumpmax, const int& M,
            Rcpp::List& jumpl, const Rcpp::List& Vparms);
void fsetjump_v(arma::vec& jump, const double& jumpmin,
                const double& jumpmax);
void fsetjump_r(arma::rowvec& jump, const double& jumpmin,
                const double& jumpmax);
void fsetjump_d(double& jump, const double& jumpmin,
                const double& jumpmax);
void cneighbor(const double& N1, const double& N2, const double& N,
               const arma::mat& trait, const arma::mat& Xnonard, 
               const int& m, arma::umat &neighbor, arma::mat& weight);

void frhononARD(arma::mat& zm, arma::vec& num, arma::vec& dm, const double& logCpzeta,  const int& N1m, const int& N2m,
                const int& Nm,  const int& P, const arma::umat& neighbor, const arma::mat& weight, 
                const arma::uvec& iARD, const arma::uvec& inonARD);
arma::mat fdnetARD(arma::mat& zm, 
                   arma::vec& num, 
                   arma::vec& dm, 
                   const int& N1m, 
                   const int& N2m, 
                   const int& Nm, 
                   const int& Pm, 
                   const double& zetam, 
                   const double& logCpzetam, 
                   const arma::umat& neighborm, 
                   const arma::mat& weightm, 
                   const arma::uvec& iARDm, 
                   const arma::uvec& inonARDm);

Rcpp::List createlistmat(const int& M, const arma::vec& Krho);
void addtolistmat(const int& M, Rcpp::List& listmat, Rcpp::List& Rho);
#endif
