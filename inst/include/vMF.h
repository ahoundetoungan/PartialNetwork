#ifndef PACKAGECPdetect_vMF_H
#define PACKAGECPdetect_vMF_H

#include <RcppArmadillo.h>

arma::mat rvMFone(const int p, const arma::vec& theta);

double logCpvMFcpp(const int& p, const double& k);

arma::mat fdnetZNU(const arma::vec& ZNU,
                   const int& Nm,
                   const int& P);
#endif
