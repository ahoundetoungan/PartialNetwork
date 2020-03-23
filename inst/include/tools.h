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
#endif
