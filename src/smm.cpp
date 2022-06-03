// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define NDEBUG 1

using namespace Rcpp;
using namespace arma;
using namespace std;

// This function draws G for one group
arma::mat fGm(const arma::mat& dnm,
              const int& Nm){
  arma::mat matunif(Nm, Nm, arma::fill::randu);
  arma::mat Gm = arma::normalise(arma::conv_to<mat>::from((matunif < dnm)), 1, 1);
  Gm.diag()    = arma::zeros(Nm);
  return Gm;
}

// Pm is the maximal power - 1 of the GMM
// This function returns Z1 as [Z1 H*Z2 H^2*Z2 ... H^Pm*Z2]
void fZ(arma::mat& Z1,
        const arma::mat& Z2,
        const arma::mat& H,
        const int& Pm){
  arma::mat tmp = Z2;
  for(int p(0); p < Pm; ++p){
    tmp = H*tmp;
    Z1  = arma::join_rows(Z1, tmp);
  }
}

// same function for fixed effects
void fZfe(arma::mat& Z1,
          const arma::mat& Z2,
          const arma::mat& H,
          const int& Pm){
  arma::mat tmp = Z2;
  for(int p(0); p < Pm; ++p){
    tmp = H*tmp;
    Z1  = arma::join_rows(Z1, tmp.each_row() - arma::mean(tmp, 0));
  }
}

// The functions starting by fgmm compute the objective function of the gmm estimator
// The functions starting by fbeta compute beta, gamma as function of alpha
// The functions starting by falbeta compute alpha, beta, gamma directly
// The functions starting by fvmzeta compute the mean and the variance of zeta
// the following number 0, 1, 2, 3  stands for Gy and GX are observed, GX is observed and not Gy, 
// Gy is observed and not GX, neither Gy nor GX is observed
// the following nc means no contextual effects
// fe means fixed effects

//Ra is R(alpha)
//Da is D(alpha)
//Day is D(alpha)y

//X1 own effects
//X2 contextual effects
//V is [X1, GX2]

//*************************** GX observed, Gy nobserved
//*** Contextual, No fixed effects
//[[Rcpp::export]]
arma::vec falbeta0(const int& R,
                   List& distr, 
                   const arma::vec& y, 
                   const arma::vec& Gy,
                   const arma::mat& GX2,
                   const arma::mat& V, 
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& Kx2, 
                   const int& ninstr,
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat Vpl =  arma::join_rows(Gy, V);
  arma::mat dG, dZ, dZVpl(ninstr, 1 + Kx1 + Kx2, arma::fill::zeros);
  arma::vec dZy(ninstr, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dG, Pm);
      dZy         += dZ.t()*ym;
      dZVpl       += dZ.t()*Vplm;
    }
  }
  arma::mat tmp    = dZVpl.t()*W;
  return arma::solve(tmp*dZVpl, tmp*dZy);
}

//[[Rcpp::export]]
List fmvzeta0(const double& alpha,
              const arma::vec& beta,
              const int& R,
              List& distr, 
              const arma::vec& y, 
              const arma::vec& Gy,
              const arma::mat& GX2,
              const arma::mat& V, 
              const arma::mat& W,
              const int& Kx1, 
              const int& Kx2, 
              const int& ninstr,
              const int& M, 
              const arma::vec& N, 
              const int& Pm,
              const arma::vec& Ncum){
  arma::mat Vpl   =  arma::join_rows(Gy, V);
  arma::vec theta = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec dZym(ninstr, arma::fill::zeros);
    arma::mat dZVplm(ninstr, 1 + Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dG, Pm);
      dZym        += dZ.t()*ym;
      dZVplm      += dZ.t()*Vplm;
    }
    matM.col(m)    = (dZym - dZVplm*theta)/R;
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH0(const double& alpha,
               const arma::vec& beta,
               const int& R,
               List& distr, 
               const arma::vec& y, 
               const arma::vec& Gy,
               const arma::mat& GX2,
               const arma::mat& V, 
               const arma::mat& W,
               const int& Kx1, 
               const int& Kx2, 
               const int& ninstr,
               const int& M, 
               const arma::vec& N, 
               const int& Pm,
               const arma::vec& Ncum){
  arma::mat Vpl   =  arma::join_rows(Gy, V);
  arma::vec theta = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ;
  arma::mat matM(ninstr, M);
  arma::mat dZVpl(ninstr, 1 + Kx1 + Kx2, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec dZym(ninstr, arma::fill::zeros);
    arma::mat dZVplm(ninstr, 1 + Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dG, Pm);
      dZym        += dZ.t()*ym;
      dZVplm      += dZ.t()*Vplm;
    }
    matM.col(m)    = (dZym - dZVplm*theta)/R;
    dZVpl         += dZVplm;
  }
  arma::mat deM    = dZVpl/(-Ncum(M)*R);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}
//*** No Contextual, No fixed effects
//[[Rcpp::export]]
arma::vec falbeta0nc(const int& R,
                     List& distr, 
                     const arma::vec& y, 
                     const arma::vec& Gy,
                     const arma::mat& GX2,
                     const arma::mat& V, 
                     const arma::mat& W,
                     const int& Kx1, 
                     const int& ninstr,
                     const int& M, 
                     const arma::vec& N, 
                     const int& Pm,
                     const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, V);
  arma::mat dG, dZ = arma::join_rows(V, GX2), dZm, dZVpl(ninstr, 1 + Kx1, arma::fill::zeros);
  arma::vec dZy(ninstr, arma::fill::zeros);
  if(Pm > 0){
    for(int m(0); m < M; ++m){
      int Nm         = N(m);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::mat dnm  = distr(m);
      arma::vec ym   = y.subvec(n1, n2);
      arma::mat GX2m = GX2.rows(n1, n2);
      arma::mat Vm   = V.rows(n1, n2);
      arma::mat Vplm = Vpl.rows(n1, n2);
      for(int r(0); r < R; ++r){
        dG           = fGm(dnm, Nm);
        arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dG, Pm);
        dZy         += dZm.t()*ym;
        dZVpl       += dZm.t()*Vplm;
      }
    }
  } else {
    dZy             = dZ.t()*y;
    dZVpl           = dZ.t()*Vpl;
  }
  
  arma::mat tmp    = dZVpl.t()*W;
  return arma::solve(tmp*dZVpl, tmp*dZy);
}

//[[Rcpp::export]]
List fmvzeta0nc(const double& alpha,
                const arma::vec& beta,
                const int& R,
                List& distr, 
                const arma::vec& y, 
                const arma::vec& Gy,
                const arma::mat& GX2,
                const arma::mat& V, 
                const arma::mat& W,
                const int& Kx1, 
                const int& ninstr,
                const int& M, 
                const arma::vec& N, 
                const int& Pm,
                const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, V);
  arma::vec theta  = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ = arma::join_rows(V, GX2), dZm;
  arma::mat matM(ninstr, M);
  if(Pm > 0){
    for(int m(0); m < M; ++m){
      arma::vec dZym(ninstr, arma::fill::zeros);
      arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
      int Nm         = N(m);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::mat dnm  = distr(m);
      arma::vec ym   = y.subvec(n1, n2);
      arma::mat GX2m = GX2.rows(n1, n2);
      arma::mat Vm   = V.rows(n1, n2);
      arma::mat Vplm = Vpl.rows(n1, n2);
      for(int r(0); r < R; ++r){
        dG           = fGm(dnm, Nm);
        arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dG, Pm);
        dZym        += dZm.t()*ym;
        dZVplm      += dZm.t()*Vplm;
      }
      matM.col(m)    = (dZym - dZVplm*theta)/R;
    }
  } else {
    for(int m(0); m < M; ++m){
      arma::vec dZym(ninstr, arma::fill::zeros);
      arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::vec ym   = y.subvec(n1, n2);
      arma::mat Vplm = Vpl.rows(n1, n2);
      arma::mat dZm  = dZ.rows(n1, n2); 
      dZym          += dZm.t()*ym;
      dZVplm        += dZm.t()*Vplm;
      matM.col(m)    = dZym - dZVplm*theta;
    }
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH0nc(const double& alpha,
                 const arma::vec& beta,
                 const int& R,
                 List& distr, 
                 const arma::vec& y, 
                 const arma::vec& Gy,
                 const arma::mat& GX2,
                 const arma::mat& V, 
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, V);
  arma::vec theta  = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ = arma::join_rows(V, GX2), dZm;
  arma::mat matM(ninstr, M);
  arma::mat dZVpl(ninstr, 1 + Kx1, arma::fill::zeros);
  if(Pm > 0){
    for(int m(0); m < M; ++m){
      arma::vec dZym(ninstr, arma::fill::zeros);
      arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
      int Nm         = N(m);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::mat dnm  = distr(m);
      arma::vec ym   = y.subvec(n1, n2);
      arma::mat GX2m = GX2.rows(n1, n2);
      arma::mat Vm   = V.rows(n1, n2);
      arma::mat Vplm = Vpl.rows(n1, n2);
      for(int r(0); r < R; ++r){
        dG           = fGm(dnm, Nm);
        arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dG, Pm);
        dZym        += dZm.t()*ym;
        dZVplm      += dZm.t()*Vplm;
      }
      matM.col(m)    = (dZym - dZVplm*theta)/R;
      dZVpl         += dZVplm;
    }
  } else {
    for(int m(0); m < M; ++m){
      arma::vec dZym(ninstr, arma::fill::zeros);
      arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::vec ym   = y.subvec(n1, n2);
      arma::mat Vplm = Vpl.rows(n1, n2);
      arma::mat dZm  = dZ.rows(n1, n2); 
      dZym          += dZm.t()*ym;
      dZVplm        += dZm.t()*Vplm;
      matM.col(m)    = dZym - dZVplm*theta;
      dZVpl         += dZVplm;
    }
  }
  arma::mat deM     = dZVpl/(-Ncum(M)*R);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//*** Contextual, fixed effects
//[[Rcpp::export]]
arma::vec falbeta0fe(const int& R,
                     List& distr, 
                     const arma::vec& y, 
                     const arma::vec& Gy,
                     const arma::mat& GX2,
                     const arma::mat& V, 
                     const arma::mat& W,
                     const int& Kx1, 
                     const int& Kx2, 
                     const int& ninstr,
                     const int& M, 
                     const arma::vec& N, 
                     const int& Pm,
                     const arma::vec& Ncum){
  arma::mat Vpl =  arma::join_rows(Gy, V);
  arma::mat dG, dZ, dZVpl(ninstr, 1 + Kx1 + Kx2, arma::fill::zeros);
  arma::vec dZy(ninstr, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dG, Pm); dZ.each_row() -= arma::mean(dZ, 0);
      dZy         += dZ.t()*ym;
      dZVpl       += dZ.t()*Vplm;
    }
  }
  arma::mat tmp    = dZVpl.t()*W;
  return arma::solve(tmp*dZVpl, tmp*dZy);
}

//[[Rcpp::export]]
List fmvzeta0fe(const double& alpha,
                const arma::vec& beta,
                const int& R,
                List& distr, 
                const arma::vec& y, 
                const arma::vec& Gy,
                const arma::mat& GX2,
                const arma::mat& V, 
                const arma::mat& W,
                const int& Kx1, 
                const int& Kx2, 
                const int& ninstr,
                const int& M, 
                const arma::vec& N, 
                const int& Pm,
                const arma::vec& Ncum){
  arma::mat Vpl   =  arma::join_rows(Gy, V);
  arma::vec theta = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec dZym(ninstr, arma::fill::zeros);
    arma::mat dZVplm(ninstr, 1 + Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dG, Pm); dZ.each_row() -= arma::mean(dZ, 0);
      dZym        += dZ.t()*ym;
      dZVplm      += dZ.t()*Vplm;
    }
    matM.col(m)    = (dZym - dZVplm*theta)/R;
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH0fe(const double& alpha,
                 const arma::vec& beta,
                 const int& R,
                 List& distr, 
                 const arma::vec& y, 
                 const arma::vec& Gy,
                 const arma::mat& GX2,
                 const arma::mat& V, 
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& Kx2, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat Vpl   =  arma::join_rows(Gy, V);
  arma::vec theta = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ;
  arma::mat matM(ninstr, M);
  arma::mat dZVpl(ninstr, 1 + Kx1 + Kx2, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec dZym(ninstr, arma::fill::zeros);
    arma::mat dZVplm(ninstr, 1 + Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dG, Pm); dZ.each_row() -= arma::mean(dZ, 0);
      dZym        += dZ.t()*ym;
      dZVplm      += dZ.t()*Vplm;
    }
    matM.col(m)    = (dZym - dZVplm*theta)/R;
    dZVpl         += dZVplm;
  }
  arma::mat deM    = dZVpl/(-Ncum(M)*R);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//*** No Contextual, fixed effects
//[[Rcpp::export]]
arma::vec falbeta0ncfe(const int& R,
                       List& distr, 
                       const arma::vec& y, 
                       const arma::vec& Gy,
                       const arma::mat& GX2,
                       const arma::mat& V, 
                       const arma::mat& W,
                       const int& Kx1, 
                       const int& ninstr,
                       const int& M, 
                       const arma::vec& N, 
                       const int& Pm,
                       const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, V);
  arma::mat dG, dZ = arma::join_rows(V, GX2), dZm, dZVpl(ninstr, 1 + Kx1, arma::fill::zeros);
  arma::vec dZy(ninstr, arma::fill::zeros);
  if(Pm > 0){
    for(int m(0); m < M; ++m){
      int Nm         = N(m);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::mat dnm  = distr(m);
      arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
      arma::mat GX2m = GX2.rows(n1, n2);
      arma::mat Vm   = V.rows(n1, n2);
      arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
      for(int r(0); r < R; ++r){
        dG           = fGm(dnm, Nm);
        arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dG, Pm); dZm.each_row() -= arma::mean(dZm, 0);
        dZy         += dZm.t()*ym;
        dZVpl       += dZm.t()*Vplm;
      }
    }
  } else {
    for(int m(0); m < M; ++m){
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
      arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
      arma::mat dZm  = dZ.rows(n1, n2); dZm.each_row() -= arma::mean(dZm, 0);
      dZy            += dZm.t()*ym;
      dZVpl          += dZm.t()*Vplm;
    }
    
  }
  arma::mat tmp    = dZVpl.t()*W;
  return arma::solve(tmp*dZVpl, tmp*dZy);
}

//[[Rcpp::export]]
List fmvzeta0ncfe(const double& alpha,
                  const arma::vec& beta,
                  const int& R,
                  List& distr, 
                  const arma::vec& y, 
                  const arma::vec& Gy,
                  const arma::mat& GX2,
                  const arma::mat& V, 
                  const arma::mat& W,
                  const int& Kx1, 
                  const int& ninstr,
                  const int& M, 
                  const arma::vec& N, 
                  const int& Pm,
                  const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, V);
  arma::vec theta  = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ = arma::join_rows(V, GX2), dZm;
  arma::mat matM(ninstr, M);
  if(Pm > 0){
    for(int m(0); m < M; ++m){
      arma::vec dZym(ninstr, arma::fill::zeros);
      arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
      int Nm         = N(m);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::mat dnm  = distr(m);
      arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
      arma::mat GX2m = GX2.rows(n1, n2);
      arma::mat Vm   = V.rows(n1, n2);
      arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
      for(int r(0); r < R; ++r){
        dG           = fGm(dnm, Nm);
        arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dG, Pm); dZm.each_row() -= arma::mean(dZm, 0);
        dZym        += dZm.t()*ym;
        dZVplm      += dZm.t()*Vplm;
      }
      matM.col(m)    = (dZym - dZVplm*theta)/R;
    }
  } else {
    for(int m(0); m < M; ++m){
      arma::vec dZym(ninstr, arma::fill::zeros);
      arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
      arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
      arma::mat dZm  = dZ.rows(n1, n2); dZm.each_row() -= arma::mean(dZm, 0);
      dZym          += dZm.t()*ym;
      dZVplm        += dZm.t()*Vplm;
      matM.col(m)    = dZym - dZVplm*theta;
    }
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH0ncfe(const double& alpha,
                   const arma::vec& beta,
                   const int& R,
                   List& distr, 
                   const arma::vec& y, 
                   const arma::vec& Gy,
                   const arma::mat& GX2,
                   const arma::mat& V, 
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& ninstr,
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, V);
  arma::vec theta  = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ = arma::join_rows(V, GX2), dZm;
  arma::mat matM(ninstr, M);
  arma::mat dZVpl(ninstr, 1 + Kx1, arma::fill::zeros);
  if(Pm > 0){
    for(int m(0); m < M; ++m){
      arma::vec dZym(ninstr, arma::fill::zeros);
      arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
      int Nm         = N(m);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::mat dnm  = distr(m);
      arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
      arma::mat GX2m = GX2.rows(n1, n2);
      arma::mat Vm   = V.rows(n1, n2);
      arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
      for(int r(0); r < R; ++r){
        dG           = fGm(dnm, Nm);
        arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dG, Pm); dZm.each_row() -= arma::mean(dZm, 0);
        dZym        += dZm.t()*ym;
        dZVplm      += dZm.t()*Vplm;
      }
      matM.col(m)    = (dZym - dZVplm*theta)/R;
      dZVpl         += dZVplm;
    }
  } else {
    for(int m(0); m < M; ++m){
      arma::vec dZym(ninstr, arma::fill::zeros);
      arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
      int n1         = Ncum(m);
      int n2         = Ncum(m + 1) - 1;
      arma::vec ym   = y.subvec(n1, n2); ym  -= mean(ym);
      arma::mat Vplm = Vpl.rows(n1, n2); Vplm.each_row() -= arma::mean(Vplm, 0);
      arma::mat dZm  = dZ.rows(n1, n2); dZm.each_row() -= arma::mean(dZm, 0);
      dZym          += dZm.t()*ym;
      dZVplm        += dZm.t()*Vplm;
      matM.col(m)    = dZym - dZVplm*theta;
      dZVpl         += dZVplm;
    }
  }
  arma::mat deM    = dZVpl/(-Ncum(M)*R);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}
//****************************************** GX observed, Gy not observed
//*** Contextual, No fixed effects
//[[Rcpp::export]]
arma::vec fbeta1(const double& alpha,
                 arma::vec& Day,
                 arma::mat& Ra,
                 const int& R,
                 const int& S,
                 const int& T, 
                 List& distr, 
                 List& Ilist, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::mat& GX2,
                 const arma::mat& V, 
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& Kx2, 
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat ddG, dA, dZ, ddGX2, ddZ, ddV, dddG;
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dddG, Pm);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm);
        Day         += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddV        = arma::join_rows(X1m, ddGX2);
          ddZ        = ddV; fZ(ddZ, ddGX2, dddG, Pm);
          Ra        += (ddZ.t()*(dA*arma::solve(Im - alpha*ddG, ddV) - ddV));
        }
        Ra          += (T*dZ.t()*Vm);
      }
    }
  }
  Day     *= T;
  return arma::solve(Ra.t()*W*Ra, Ra.t()*(W*Day));
}
//[[Rcpp::export]]
double fgmm1(const double& alpha,
             const int& R,
             const int& S,
             const int& T, 
             List& distr, 
             List& Ilist, 
             const arma::vec& y, 
             const arma::mat& X1, 
             const arma::mat& X2, 
             const arma::mat& GX2,
             const arma::mat& V, 
             const arma::mat& W,
             const int& Kx1, 
             const int& Kx2, 
             const int& ninstr,
             const int& M, 
             const arma::vec& N, 
             const int& Pm,
             const arma::vec& Ncum){
  arma::vec Day(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::vec beta = fbeta1(alpha, Day, Ra, R, S, T, distr, Ilist, y, X1, X2, GX2, V, W, Kx1, Kx2, M, N, Pm, Ncum);
  arma::vec h = Day - Ra*beta;
  return arma::dot(h, W*h)/(Ncum(M)*Ncum(M)*R*R*T*T*S*S);
}

//[[Rcpp::export]]
List fmvzeta1(const double& alpha,
              const arma::vec& beta,
              const int& R,
              const int& S,
              const int& T, 
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::mat& X1, 
              const arma::mat& X2, 
              const arma::mat& GX2,
              const arma::mat& V, 
              const arma::mat& W,
              const int& Kx1, 
              const int& Kx2, 
              const int& ninstr,
              const int& M, 
              const arma::vec& N, 
              const int& Pm,
              const arma::vec& Ncum){
  arma::mat ddG, dA, dZ, ddGX2, ddZ, ddV, dddG;
  arma::mat matM(ninstr, M);
  
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dddG, Pm);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm);
        Daym        += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddV        = arma::join_rows(X1m, ddGX2);
          ddZ        = ddV; fZ(ddZ, ddGX2, dddG, Pm);
          Ram       += (ddZ.t()*(dA*arma::solve(Im - alpha*ddG, ddV) - ddV));
        }
        Ram         += (T*dZ.t()*Vm);
      }
    }
    matM.col(m)    = Daym/(R*S) - Ram*beta/(R*S*T);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH1(const double& alpha,
               const arma::vec& beta,
               const int& R,
               const int& S,
               const int& T, 
               List& distr, 
               List& Ilist, 
               const arma::vec& y, 
               const arma::mat& X1, 
               const arma::mat& X2, 
               const arma::mat& GX2,
               const arma::mat& V, 
               const arma::mat& W,
               const int& Kx1, 
               const int& Kx2, 
               const int& ninstr,
               const int& M, 
               const arma::vec& N, 
               const int& Pm,
               const arma::vec& Ncum){
  arma::mat dG, ddG, dA, dZ, ddGX2, ddZ, ddA, ddV, ddAV, dddG;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1 + Kx2, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dddG, Pm);
      for(int s(0); s < S; ++s){
        dG           = fGm(dnm, Nm);
        dA           = Im - alpha*dG;
        Daym        += (dZ.t()*dA*ym);
        Dgym        += (dZ.t()*dG*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddV        = arma::join_rows(X1m, ddGX2);
          ddZ        = ddV; fZ(ddZ, ddGX2, dddG, Pm);
          ddA        = Im - alpha*ddG;
          ddAV       = arma::solve(ddA, ddV);
          Ram       += (ddZ.t()*(dA*ddAV - ddV));
          tmp1      += (ddZ.t()*dG*ddAV);
          tmp2      += (ddZ.t()*dA*arma::solve(ddA, ddG*ddAV));
        }
        Ram         += (T*dZ.t()*Vm);
      }
    }
    matM.col(m)    = Daym/(R*S) - (Ram*beta)/(R*S*T);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + Kx2 + 1);
  deM.col(0)             = (T*Dgy - (tmp1 - tmp2)*beta)/(-Ncum(M)*R*S*T);
  deM.cols(1, Kx1 + Kx2) = Ra/(-Ncum(M)*R*S*T);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//*** No Contextual, No fixed effects
//[[Rcpp::export]]
arma::vec fbeta1nc(const double& alpha,
                   arma::vec& Day,
                   arma::mat& Ra,
                   const int& R,
                   const int& S,
                   const int& T, 
                   List& distr, 
                   List& Ilist, 
                   const arma::vec& y, 
                   const arma::mat& X1, 
                   const arma::mat& X2, 
                   const arma::mat& GX2,
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ = arma::join_rows(X1, GX2), ddGX2, ddZ;
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dddG, Pm);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm);
        Day         += (dZm.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddZ        = arma::join_rows(X1m, ddGX2); fZ(ddZ, ddGX2, dddG, Pm);
          Ra        += (ddZ.t()*(dA*arma::solve(Im - alpha*ddG, X1m) - X1m));
        }
        Ra          += (T*dZm.t()*X1m);
      }
    }
  }
  Day     *= T;
  return arma::solve(Ra.t()*W*Ra, Ra.t()*(W*Day));
}


//[[Rcpp::export]]
double fgmm1nc(const double& alpha,
               const int& R,
               const int& S,
               const int& T, 
               List& distr, 
               List& Ilist, 
               const arma::vec& y, 
               const arma::mat& X1, 
               const arma::mat& X2, 
               const arma::mat& GX2,
               const arma::mat& W,
               const int& Kx1, 
               const int& ninstr,
               const int& M, 
               const arma::vec& N, 
               const int& Pm,
               const arma::vec& Ncum){
  arma::vec Day(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1, arma::fill::zeros);
  arma::vec beta = fbeta1nc(alpha, Day, Ra, R, S, T, distr, Ilist, y, X1, X2, GX2, W, Kx1, M, N, Pm, Ncum);
  arma::vec h = Day - Ra*beta;
  return arma::dot(h, W*h)/(Ncum(M)*Ncum(M)*R*R*T*T*S*S);
}

//[[Rcpp::export]]
List fmvzeta1nc(const double& alpha,
                const arma::vec& beta,
                const int& R,
                const int& S,
                const int& T, 
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::mat& X1, 
                const arma::mat& X2, 
                const arma::mat& GX2,
                const arma::mat& W,
                const int& Kx1, 
                const int& ninstr,
                const int& M, 
                const arma::vec& N, 
                const int& Pm,
                const arma::vec& Ncum){
  arma::vec X1b = X1*beta.head(Kx1);
  arma::mat dddG, ddG, dA, dZ = arma::join_rows(X1, GX2), ddGX2, ddZ;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Rabm(ninstr, Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec X1bm = X1b.subvec(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dddG, Pm);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm);
        Daym        += (dZm.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddZ        = arma::join_rows(X1m, ddGX2); fZ(ddZ, ddGX2, dddG, Pm);
          Rabm      += (ddZ.t()*(dA*arma::solve(Im - alpha*ddG, X1bm) - X1bm));
        }
        Rabm        += (T*dZm.t()*X1bm);
      }
    }
    matM.col(m)    = Daym/(R*S) - Rabm/(R*S*T);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH1nc(const double& alpha,
                 const arma::vec& beta,
                 const int& R,
                 const int& S,
                 const int& T, 
                 List& distr, 
                 List& Ilist, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::mat& GX2,
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat dG, ddG, dA, dZ = arma::join_rows(X1, GX2), ddGX2, ddZ, ddA, ddV, ddAV, dddG;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dddG, Pm);
      for(int s(0); s < S; ++s){
        dG           = fGm(dnm, Nm);
        dA           = Im - alpha*dG;
        Daym        += (dZm.t()*dA*ym);
        Dgym        += (dZm.t()*dG*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddZ        = arma::join_rows(X1m, ddGX2); fZ(ddZ, ddGX2, dddG, Pm);
          ddA        = Im - alpha*ddG;
          ddAV       = arma::solve(ddA, X1m);
          Ram       += (ddZ.t()*(dA*ddAV - X1m));
          tmp1      += (ddZ.t()*dG*ddAV);
          tmp2      += (ddZ.t()*dA*arma::solve(ddA, ddG*ddAV));
        }
        Ram         += (T*dZm.t()*X1m);
      }
    }
    matM.col(m)    = Daym/(R*S) - (Ram*beta)/(R*S*T);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + 1);
  deM.col(0)       = (T*Dgy - (tmp1 - tmp2)*beta)/(-Ncum(M)*R*S*T);
  deM.cols(1, Kx1) = Ra/(-Ncum(M)*R*S*T);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}
//*** Contextual, fixed effects
//[[Rcpp::export]]
arma::vec fbeta1fe(const double& alpha,
                   arma::vec& Day,
                   arma::mat& Ra,
                   const int& R,
                   const int& S,
                   const int& T, 
                   List& distr, 
                   List& Ilist, 
                   const arma::vec& y, 
                   const arma::mat& X1, 
                   const arma::mat& X2, 
                   const arma::mat& GX2,
                   const arma::mat& V, 
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& Kx2, 
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ, ddGX2, ddZ, ddV;
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2); 
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2); 
    arma::mat Vmc  = Vm.each_row() - mean(Vm, 0); 
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dddG, Pm); dZ.each_row() -= arma::mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm); dA.each_row() -= arma::mean(dA, 0);
        Day         += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddV        = arma::join_rows(X1m, ddGX2); 
          ddZ        = ddV; fZ(ddZ, ddGX2, dddG, Pm); ddZ.each_row() -= arma::mean(ddZ, 0);
          Ra        += (ddZ.t()*(dA*arma::solve(Im - alpha*ddG, ddV) - (ddV.each_row() - mean(ddV, 0))));
        }
        Ra          += (T*dZ.t()*Vmc);
      }
    }
  }
  Day     *= T;
  return arma::solve(Ra.t()*W*Ra, Ra.t()*(W*Day));
}

//[[Rcpp::export]]
double fgmm1fe(const double& alpha,
               const int& R,
               const int& S,
               const int& T, 
               List& distr, 
               List& Ilist, 
               const arma::vec& y, 
               const arma::mat& X1, 
               const arma::mat& X2, 
               const arma::mat& GX2,
               const arma::mat& V, 
               const arma::mat& W,
               const int& Kx1, 
               const int& Kx2, 
               const int& ninstr,
               const int& M, 
               const arma::vec& N, 
               const int& Pm,
               const arma::vec& Ncum){
  arma::vec Day(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::vec beta = fbeta1fe(alpha, Day, Ra, R, S, T, distr, Ilist, y, X1, X2, GX2, V, W, Kx1, Kx2, M, N, Pm, Ncum);
  arma::vec h = Day - Ra*beta;
  return arma::dot(h, W*h)/(Ncum(M)*Ncum(M)*R*R*T*T*S*S);
}

//[[Rcpp::export]]
List fmvzeta1fe(const double& alpha,
                const arma::vec& beta,
                const int& R,
                const int& S,
                const int& T, 
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::mat& X1, 
                const arma::mat& X2, 
                const arma::mat& GX2,
                const arma::mat& V, 
                const arma::mat& W,
                const int& Kx1, 
                const int& Kx2, 
                const int& ninstr,
                const int& M, 
                const arma::vec& N, 
                const int& Pm,
                const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ, ddGX2, ddZ, ddV;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2); 
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2); 
    arma::mat Vmc  = Vm.each_row() - mean(Vm, 0); 
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dddG, Pm); dZ.each_row() -= arma::mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm); dA.each_row() -= arma::mean(dA, 0);
        Daym        += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddV        = arma::join_rows(X1m, ddGX2); 
          ddZ        = ddV; fZ(ddZ, ddGX2, dddG, Pm); ddZ.each_row() -= arma::mean(ddZ, 0);
          Ram       += (ddZ.t()*(dA*arma::solve(Im - alpha*ddG, ddV) - (ddV.each_row() - mean(ddV, 0))));
        }
        Ram         += (T*dZ.t()*Vmc);
      }
    }
    matM.col(m)    = Daym/(R*S) - Ram*beta/(R*S*T);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH1fe(const double& alpha,
                 const arma::vec& beta,
                 const int& R,
                 const int& S,
                 const int& T, 
                 List& distr, 
                 List& Ilist, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::mat& GX2,
                 const arma::mat& V, 
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& Kx2, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat dG, ddG, dA, dZ, ddGX2, ddZ, ddV, ddA, ddAV, dGc, dddG;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1 + Kx2, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2); 
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    arma::mat Vm   = V.rows(n1, n2); 
    arma::mat Vmc  = Vm.each_row() - mean(Vm, 0); 
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dZ           = Vm;  fZ(dZ, GX2m, dddG, Pm); dZ.each_row() -= arma::mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dG           = fGm(dnm, Nm);
        dGc          = dG.each_row() - mean(dG, 0);
        dA           = Im - alpha*dG; dA.each_row() -= arma::mean(dA, 0);
        Daym        += (dZ.t()*dA*ym);
        Dgym        += (dZ.t()*dGc*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddV        = arma::join_rows(X1m, ddGX2); 
          ddZ        = ddV; fZ(ddZ, ddGX2, dddG, Pm); ddZ.each_row() -= arma::mean(ddZ, 0);
          ddA        = Im - alpha*ddG;
          ddAV       = arma::solve(ddA, ddV);
          Ram       += (ddZ.t()*(dA*ddAV - (ddV.each_row() - mean(ddV, 0))));
          tmp1      += (ddZ.t()*dGc*ddAV);
          tmp2      += (ddZ.t()*dA*arma::solve(ddA, ddG*ddAV));
        }
        Ram         += (T*dZ.t()*Vmc);
      }
    }
    matM.col(m)    = Daym/(R*S) - (Ram*beta)/(R*S*T);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + Kx2 + 1);
  deM.col(0)             = (T*Dgy - (tmp1 - tmp2)*beta)/(-Ncum(M)*R*S*T);
  deM.cols(1, Kx1 + Kx2) = Ra/(-Ncum(M)*S*T);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//*** No Contextual, fixed effects
//[[Rcpp::export]]
arma::vec fbeta1ncfe(const double& alpha,
                     arma::vec& Day,
                     arma::mat& Ra,
                     const int& R,
                     const int& S,
                     const int& T, 
                     List& distr, 
                     List& Ilist, 
                     const arma::vec& y, 
                     const arma::mat& X1, 
                     const arma::mat& X2, 
                     const arma::mat& GX2,
                     const arma::mat& W,
                     const int& Kx1, 
                     const int& M, 
                     const arma::vec& N, 
                     const int& Pm,
                     const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ = arma::join_rows(X1, GX2), ddGX2, ddZ;
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X1mc = X1m.each_row() - arma::mean(X1m, 0);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dddG, Pm); dZm.each_row() -= arma::mean(dZm, 0);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm); dA.each_row() -= arma::mean(dA, 0);
        Day         += (dZm.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddZ        = arma::join_rows(X1m, ddGX2); fZ(ddZ, ddGX2, dddG, Pm); ddZ.each_row() -= arma::mean(ddZ, 0);
          Ra        += (ddZ.t()*(dA*arma::solve(Im - alpha*ddG, X1m) - X1mc));
        }
        Ra          += (T*dZm.t()*X1mc);
      }
    }
  }
  Day     *= T;
  return arma::solve(Ra.t()*W*Ra, Ra.t()*(W*Day));
}

//[[Rcpp::export]]
double fgmm1ncfe(const double& alpha,
                 const int& R,
                 const int& S,
                 const int& T, 
                 List& distr, 
                 List& Ilist, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::mat& GX2,
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::vec Day(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1, arma::fill::zeros);
  arma::vec beta = fbeta1ncfe(alpha, Day, Ra, R, S, T, distr, Ilist, y, X1, X2, GX2, W, Kx1, M, N, Pm, Ncum);
  arma::vec h = Day - Ra*beta;
  return arma::dot(h, W*h)/(Ncum(M)*Ncum(M)*R*R*T*T*S*S);
}

//[[Rcpp::export]]
List fmvzeta1ncfe(const double& alpha,
                  const arma::vec& beta,
                  const int& R,
                  const int& S,
                  const int& T, 
                  List& distr, 
                  List& Ilist, 
                  const arma::vec& y, 
                  const arma::mat& X1, 
                  const arma::mat& X2, 
                  const arma::mat& GX2,
                  const arma::mat& W,
                  const int& Kx1, 
                  const int& ninstr,
                  const int& M, 
                  const arma::vec& N, 
                  const int& Pm,
                  const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ = arma::join_rows(X1, GX2), ddGX2, ddZ;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X1mc = X1m.each_row() - arma::mean(X1m, 0);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      for(int s(0); s < S; ++s){
        dddG         = fGm(dnm, Nm);
        dA           = Im - alpha*fGm(dnm, Nm); dA.each_row() -= arma::mean(dA, 0);
        arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dddG, Pm); dZm.each_row() -= arma::mean(dZm, 0);
        Daym        += (dZm.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddZ        = arma::join_rows(X1m, ddGX2); fZ(ddZ, ddGX2, dddG, Pm); ddZ.each_row() -= arma::mean(ddZ, 0);
          Ram       += (ddZ.t()*(dA*arma::solve(Im - alpha*ddG, X1m) - X1mc));
        }
        Ram         += (T*dZm.t()*X1mc);
      }
    }
    matM.col(m)    = Daym/(R*S) - Ram*beta/(R*S*T);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH1ncfe(const double& alpha,
                   const arma::vec& beta,
                   const int& R,
                   const int& S,
                   const int& T, 
                   List& distr, 
                   List& Ilist, 
                   const arma::vec& y, 
                   const arma::mat& X1, 
                   const arma::mat& X2, 
                   const arma::mat& GX2,
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& ninstr,
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat dG, ddG, dA, dZ = arma::join_rows(X1, GX2), ddGX2, ddZ, ddA, ddV, ddAV, dGc, dddG;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X1mc = X1m.each_row() - arma::mean(X1m, 0);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat GX2m = GX2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      arma::mat dZm= dZ.rows(n1, n2); fZ(dZm, GX2m, dddG, Pm); dZm.each_row() -= arma::mean(dZm, 0);
      for(int s(0); s < S; ++s){
        dG           = fGm(dnm, Nm);
        dGc          = dG.each_row() - mean(dG, 0);
        dA           = Im - alpha*dG; dA.each_row() -= arma::mean(dA, 0);
        Daym        += (dZm.t()*dA*ym);
        Dgym        += (dZm.t()*dGc*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddGX2      = ddG*X2m;
          ddZ        = arma::join_rows(X1m, ddGX2); fZ(ddZ, ddGX2, dddG, Pm); ddZ.each_row() -= arma::mean(ddZ, 0);
          ddA        = Im - alpha*ddG;
          ddAV       = arma::solve(ddA, X1m);
          Ram       += (ddZ.t()*(dA*ddAV - X1mc));
          tmp1      += (ddZ.t()*dGc*ddAV);
          tmp2      += (ddZ.t()*dA*arma::solve(ddA, ddG*ddAV));
        }
        Ram         += (T*dZm.t()*X1mc);
      }
    }
    matM.col(m)    = Daym/(R*S) - (Ram*beta)/(R*S*T);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + 1);
  deM.col(0)       = (T*Dgy - (tmp1 - tmp2)*beta)/(-Ncum(M)*R*S*T);
  deM.cols(1, Kx1) = Ra/(-Ncum(M)*R*S*T);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//****************************************** GX not observed, Gy observed
//*** Contextual, No fixed effects
//[[Rcpp::export]]
arma::vec falbeta2(const int& R,
                 const int& S, 
                 List& distr, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::vec& Gy,
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& Kx2, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat dG, dZ, dGX2, dZVpl(ninstr, Kx1 + Kx2 + 1, arma::fill::zeros);
  arma::vec dZy(ninstr, arma::fill::zeros);
  arma::mat X1p = arma::join_rows(Gy, X1);
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X1pm = X1p.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      dZ           = arma::join_rows(X1m, dGX2);  fZ(dZ, dGX2, dG, Pm);
      dZy         += (dZ.t()*ym);
      for(int s(0); s < S; ++s){
        dZVpl     += (dZ.t()*arma::join_rows(X1pm, fGm(dnm, Nm)*X2m));
      }
    }
  }
  arma::mat tmp    = dZVpl.t()*W;
  return arma::solve(tmp*dZVpl, S*tmp*dZy);
}


//[[Rcpp::export]]
List fmvzeta2(const double& alpha,
              const arma::vec& beta,
              const int& R,
              const int& S, 
              List& distr, 
              const arma::vec& y, 
              const arma::mat& X1, 
              const arma::mat& X2, 
              const arma::vec& Gy,
              const arma::mat& W,
              const int& Kx1, 
              const int& Kx2, 
              const int& ninstr,
              const int& M, 
              const arma::vec& N, 
              const int& Pm,
              const arma::vec& Ncum){
  arma::mat dG, ddG, dZ, dGX2;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec Gym  = Gy.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      dZ           = arma::join_rows(X1m, dGX2);  fZ(dZ, dGX2, dG, Pm);
      Daym        += (dZ.t()*(ym - alpha*Gym));
      for(int s(0); s < S; ++s){
        ddG        = fGm(dnm, Nm);
        Ram       += (dZ.t()*arma::join_rows(X1m, ddG*X2m));
      }
    }
    matM.col(m)    = Daym/R - Ram*beta/(R*S);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH2(const double& alpha,
               const arma::vec& beta,
               const int& R,
               const int& S, 
               List& distr, 
               const arma::vec& y, 
               const arma::mat& X1, 
               const arma::mat& X2, 
               const arma::vec& Gy,
               const arma::mat& W,
               const int& Kx1, 
               const int& Kx2, 
               const int& ninstr,
               const int& M, 
               const arma::vec& N, 
               const int& Pm,
               const arma::vec& Ncum){
  arma::mat dG, ddG, dZ, dGX2, ddV;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec Gym  = Gy.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      dZ           = arma::join_rows(X1m, dGX2);  fZ(dZ, dGX2, dG, Pm);
      Daym        += (dZ.t()*(ym - alpha*Gym));
      Dgym        += (dZ.t()*dG*ym);
      for(int s(0); s < S; ++s){
        ddG        = fGm(dnm, Nm);
        ddV        = arma::join_rows(X1m, ddG*X2m);
        Ram       += (dZ.t()*ddV);
      }
    }
    matM.col(m)    = Daym/R - Ram*beta/(R*S);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + Kx2 + 1);
  deM.col(0)             = Dgy/(-Ncum(M)*R);
  deM.cols(1, Kx1 + Kx2) = Ra/(-Ncum(M)*R*S);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}


//*** No Contextual, No fixed effects
//[[Rcpp::export]]
arma::vec falbeta2nc(const int& R,
                     List& distr, 
                     const arma::vec& y, 
                     const arma::mat& X1, 
                     const arma::mat& X2, 
                     const arma::vec& Gy,
                     const arma::mat& W,
                     const int& Kx1, 
                     const int& ninstr,
                     const int& M, 
                     const arma::vec& N, 
                     const int& Pm,
                     const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, X1);
  arma::mat dG, dZ, dGX2, dZVpl(ninstr, 1 + Kx1, arma::fill::zeros);
  arma::vec dZy(ninstr, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      arma::mat dZ = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm);
      dZy         += dZ.t()*ym;
      dZVpl       += dZ.t()*Vplm;
    }
  }
  arma::mat tmp    = dZVpl.t()*W;
  return arma::solve(tmp*dZVpl, tmp*dZy);
}

//[[Rcpp::export]]
List fmvzeta2nc(const double& alpha,
                const arma::vec& beta,
                const int& R,
                List& distr, 
                const arma::vec& y, 
                const arma::mat& X1, 
                const arma::mat& X2, 
                const arma::vec& Gy,
                const arma::mat& W,
                const int& Kx1, 
                const int& ninstr,
                const int& M, 
                const arma::vec& N, 
                const int& Pm,
                const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, X1);
  arma::vec theta  = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ, dGX2;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec dZym(ninstr, arma::fill::zeros);
    arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      arma::mat dZ = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm);
      dZym        += dZ.t()*ym;
      dZVplm      += dZ.t()*Vplm;
    }
    matM.col(m)    = (dZym - dZVplm*theta)/R;
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH2nc(const double& alpha,
                 const arma::vec& beta,
                 const int& R,
                 List& distr, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::vec& Gy,
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, X1);
  arma::vec theta  = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ, dGX2;
  arma::mat matM(ninstr, M);
  arma::mat dZVpl(ninstr, 1 + Kx1, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec dZym(ninstr, arma::fill::zeros);
    arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      arma::mat dZ = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm);
      dZym        += dZ.t()*ym;
      dZVplm      += dZ.t()*Vplm;
    }
    matM.col(m)    = (dZym - dZVplm*theta)/R;
    dZVpl         += dZVplm;
  }
  arma::mat deM    = dZVpl/(-Ncum(M)*R);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//*** Contextual, fixed effects
//[[Rcpp::export]]
arma::vec falbeta2fe(const int& R,
                   const int& S, 
                   List& distr, 
                   const arma::vec& y, 
                   const arma::mat& X1, 
                   const arma::mat& X2, 
                   const arma::vec& Gy,
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& Kx2, 
                   const int& ninstr,
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat dG, dZ, dGX2, Vpl, dZVpl(ninstr, Kx1 + Kx2 + 1, arma::fill::zeros);
  arma::vec dZy(ninstr, arma::fill::zeros);
  arma::mat X1p = arma::join_rows(Gy, X1);
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2); ym -= mean(ym);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X1pm = X1p.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      dZ           = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm); dZ.each_row() -= mean(dZ, 0);
      dZy         += (dZ.t()*ym);
      for(int s(0); s < S; ++s){
        Vpl        = arma::join_rows(X1pm, fGm(dnm, Nm)*X2m); Vpl.each_row() -= mean(Vpl, 0);
        dZVpl     += (dZ.t()*Vpl);
      }
    }
  }

  arma::mat tmp    = dZVpl.t()*W;
  return arma::solve(tmp*dZVpl, S*tmp*dZy);
}

//[[Rcpp::export]]
List fmvzeta2fe(const double& alpha,
                const arma::vec& beta,
                const int& R,
                const int& S, 
                List& distr, 
                const arma::vec& y, 
                const arma::mat& X1, 
                const arma::mat& X2, 
                const arma::vec& Gy,
                const arma::mat& W,
                const int& Kx1, 
                const int& Kx2, 
                const int& ninstr,
                const int& M, 
                const arma::vec& N, 
                const int& Pm,
                const arma::vec& Ncum){
  arma::mat dG, ddG, dZ, dGX2, tmp;
  arma::vec Ay;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec Gym  = Gy.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      dZ           = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm); dZ.each_row() -= mean(dZ, 0);
      Ay           = ym - alpha*Gym; Ay -= mean(Ay);
      Daym        += (dZ.t()*Ay);
      for(int s(0); s < S; ++s){
        ddG        = fGm(dnm, Nm);
        tmp        = arma::join_rows(X1m, ddG*X2m); tmp.each_row() -= mean(tmp, 0);
        Ram       += (dZ.t()*tmp);
      }
    }
    matM.col(m)    = Daym/R - Ram*beta/(R*S);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}


//[[Rcpp::export]]
List fmvzetaH2fe(const double& alpha,
                 const arma::vec& beta,
                 const int& R,
                 const int& S, 
                 List& distr, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::vec& Gy,
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& Kx2, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat dG, ddG, dZ, dGX2, tmp, dGc, ddGX2;
  arma::vec Ay;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1 + Kx2, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec Gym  = Gy.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGc          = dG.each_row() - mean(dG, 0);
      dGX2         = dG*X2m;
      dZ           = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm); dZ.each_row() -= mean(dZ, 0);
      Ay           = ym - alpha*Gym; Ay -= mean(Ay);
      Daym        += (dZ.t()*Ay);
      Dgym        += (dZ.t()*dGc*ym);
      for(int s(0); s < S; ++s){
        ddG        = fGm(dnm, Nm);
        ddGX2      = ddG*X2m;
        tmp        = arma::join_rows(X1m, ddGX2); tmp.each_row() -= mean(tmp, 0);
        Ram       += (dZ.t()*tmp);
      }
    }
    matM.col(m)    = Daym/R - Ram*beta/(R*S);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + Kx2 + 1);
  deM.col(0)             = Dgy/(-Ncum(M)*R);
  deM.cols(1, Kx1 + Kx2) = Ra/(-Ncum(M)*R*S);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//*** No Contextual, fixed effects
//[[Rcpp::export]]
arma::vec falbeta2ncfe(const int& R,
                       List& distr, 
                       const arma::vec& y, 
                       const arma::mat& X1, 
                       const arma::mat& X2, 
                       const arma::vec& Gy,
                       const arma::mat& W,
                       const int& Kx1, 
                       const int& ninstr,
                       const int& M, 
                       const arma::vec& N, 
                       const int& Pm,
                       const arma::vec& Ncum){
  arma::mat Vpl    =  arma::join_rows(Gy, X1);
  arma::mat dG, dZ, dGX2, dZVpl(ninstr, 1 + Kx1, arma::fill::zeros);
  arma::vec dZy(ninstr, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2); ym -= mean(ym);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2); Vpl.each_row() -= mean(Vpl, 0);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      arma::mat dZ = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm); dZ.each_row() -= mean(dZ, 0);
      dZy         += dZ.t()*ym;
      dZVpl       += dZ.t()*Vplm;
    }
  }
  arma::mat tmp    = dZVpl.t()*W;
  return arma::solve(tmp*dZVpl, tmp*dZy);
}


//[[Rcpp::export]]
List fmvzeta2ncfe(const double& alpha,
                  const arma::vec& beta,
                  const int& R,
                  List& distr, 
                  const arma::vec& y, 
                  const arma::mat& X1, 
                  const arma::mat& X2, 
                  const arma::vec& Gy,
                  const arma::mat& W,
                  const int& Kx1, 
                  const int& ninstr,
                  const int& M, 
                  const arma::vec& N, 
                  const int& Pm,
                  const arma::vec& Ncum){
  arma::mat Vpl    = arma::join_rows(Gy, X1);
  arma::vec theta  = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ, dGX2;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec dZym(ninstr, arma::fill::zeros);
    arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2); ym -= mean(ym);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2); Vpl.each_row() -= mean(Vpl, 0);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      arma::mat dZ = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm); dZ.each_row() -= mean(dZ, 0);
      dZym        += dZ.t()*ym;
      dZVplm      += dZ.t()*Vplm;
    }
    matM.col(m)    = (dZym - dZVplm*theta)/R;
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH2ncfe(const double& alpha,
                   const arma::vec& beta,
                   const int& R,
                   List& distr, 
                   const arma::vec& y, 
                   const arma::mat& X1, 
                   const arma::mat& X2, 
                   const arma::vec& Gy,
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& ninstr,
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat Vpl    = arma::join_rows(Gy, X1);
  arma::vec theta  = arma::join_cols(alpha*arma::ones(1), beta);
  arma::mat dG, dZ, dGX2;
  arma::mat matM(ninstr, M);
  arma::mat dZVpl(ninstr, 1 + Kx1, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec dZym(ninstr, arma::fill::zeros);
    arma::mat dZVplm(ninstr, 1 + Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::vec ym   = y.subvec(n1, n2); ym -= mean(ym);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::mat Vplm = Vpl.rows(n1, n2); Vpl.each_row() -= mean(Vpl, 0);
    for(int r(0); r < R; ++r){
      dG           = fGm(dnm, Nm);
      dGX2         = dG*X2m;
      arma::mat dZ = arma::join_rows(X1m, dGX2); fZ(dZ, dGX2, dG, Pm); dZ.each_row() -= mean(dZ, 0);
      dZym        += dZ.t()*ym;
      dZVplm      += dZ.t()*Vplm;
    }
    matM.col(m)    = (dZym - dZVplm*theta)/R;
    dZVpl         += dZVplm;
  }
  arma::mat deM    = dZVpl/(-Ncum(M)*R);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}
//***************************************** GX not observed, Gy not observed
//*** Contextual, No fixed effects
//[[Rcpp::export]]
arma::vec fbeta3(const double& alpha,
                 arma::vec& Day,
                 arma::mat& Ra,
                 const int& R,
                 const int& S,
                 const int& T, 
                 List& distr, 
                 List& Ilist, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& Kx2, 
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ, dddGX2;
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dddGX2       = dddG*X2m;
      dZ           = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm);
        Day         += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          Ra        += (dZ.t()*dA*arma::solve(Im - alpha*ddG, arma::join_rows(X1m, ddG*X2m)));
        }
      }
    }
  }
  Day     *= T;
  return arma::solve(Ra.t()*W*Ra, Ra.t()*(W*Day));
}

//[[Rcpp::export]]
double fgmm3(const double& alpha,
             const int& R,
             const int& S,
             const int& T, 
             List& distr, 
             List& Ilist, 
             const arma::vec& y, 
             const arma::mat& X1, 
             const arma::mat& X2, 
             const arma::mat& W,
             const int& Kx1, 
             const int& Kx2, 
             const int& ninstr,
             const int& M, 
             const arma::vec& N, 
             const int& Pm,
             const arma::vec& Ncum){
  arma::vec Day(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::vec beta = fbeta3(alpha, Day, Ra, R, S, T, distr, Ilist, y, X1, X2, W, Kx1, Kx2, M, N, Pm, Ncum);
  arma::vec h = Day - Ra*beta;
  return arma::dot(h, W*h)/(Ncum(M)*Ncum(M)*R*R*T*T*S*S);
}

//[[Rcpp::export]]
List fmvzeta3(const double& alpha,
              const arma::vec& beta,
              const int& R,
              const int& S,
              const int& T, 
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::mat& X1, 
              const arma::mat& X2, 
              const arma::mat& W,
              const int& Kx1, 
              const int& Kx2, 
              const int& ninstr,
              const int& M, 
              const arma::vec& N, 
              const int& Pm,
              const arma::vec& Ncum){
  arma::vec X1b = X1*beta.head(Kx1);
  arma::vec X2b = X2*beta.tail(Kx2);
  arma::mat dddG, ddG, dA, dZ, dddGX2;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Rabm(ninstr, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec X1bm = X1b.subvec(n1, n2);
    arma::vec X2bm = X2b.subvec(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dddGX2       = dddG*X2m;
      dZ           = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm);
        Daym        += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          Rabm      += (dZ.t()*dA*arma::solve(Im - alpha*ddG, X1bm + ddG*X2bm));
        }
      }
    }
    matM.col(m)    = Daym/(R*S) - Rabm/(R*S*T);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH3(const double& alpha,
               const arma::vec& beta,
               const int& R,
               const int& S,
               const int& T, 
               List& distr, 
               List& Ilist, 
               const arma::vec& y, 
               const arma::mat& X1, 
               const arma::mat& X2, 
               const arma::mat& W,
               const int& Kx1, 
               const int& Kx2, 
               const int& ninstr,
               const int& M, 
               const arma::vec& N, 
               const int& Pm,
               const arma::vec& Ncum){
  arma::mat dG, ddG, dddG, dA, dZ, dddGX2, ddA, ddAV;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1 + Kx2, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG         = fGm(dnm, Nm);
      dddGX2       = dddG*X2m;
      dZ           = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm);
      for(int s(0); s < S; ++s){
        dG           = fGm(dnm, Nm);
        dA           = Im - alpha*dG;
        Daym        += (dZ.t()*dA*ym);
        Dgym        += (dZ.t()*dG*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddA        = Im - alpha*ddG;
          ddAV       = arma::solve(ddA, arma::join_rows(X1m, ddG*X2m));
          tmp1      += (dZ.t()*dG*ddAV);
          tmp2      += (dZ.t()*dA*arma::solve(ddA, ddG*ddAV));
          Ram       += (dZ.t()*dA*ddAV);
        }
      }
    }
    matM.col(m)    = Daym/(R*S) - (Ram*beta)/(R*S*T);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + Kx2 + 1);
  deM.col(0)             = (T*Dgy - (tmp1 - tmp2)*beta)/(-Ncum(M)*R*S*T);
  deM.cols(1, Kx1 + Kx2) = Ra/(-Ncum(M)*R*S*T);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//*** No Contextual, No fixed effects
//[[Rcpp::export]]
arma::vec fbeta3nc(const double& alpha,
                   arma::vec& Day,
                   arma::mat& Ra,
                   const int& R,
                   const int& S,
                   const int& T, 
                   List& distr, 
                   List& Ilist, 
                   const arma::vec& y, 
                   const arma::mat& X1, 
                   const arma::mat& X2, 
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ, dddGX2;
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm);
        Day         += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          Ra        += (dZ.t()*dA*arma::solve(Im - alpha*ddG, X1m));
        }
      }
    }
  }
  Day     *= T;
  return arma::solve(Ra.t()*W*Ra, Ra.t()*(W*Day));
}

//[[Rcpp::export]]
double fgmm3nc(const double& alpha,
               const int& R,
               const int& S,
               const int& T, 
               List& distr, 
               List& Ilist, 
               const arma::vec& y, 
               const arma::mat& X1, 
               const arma::mat& X2, 
               const arma::mat& W,
               const int& Kx1, 
               const int& ninstr,
               const int& M, 
               const arma::vec& N, 
               const int& Pm,
               const arma::vec& Ncum){
  arma::vec Day(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1, arma::fill::zeros);
  arma::vec beta = fbeta3nc(alpha, Day, Ra, R, S, T, distr, Ilist, y, X1, X2, W, Kx1, M, N, Pm, Ncum);
  arma::vec h = Day - Ra*beta;
  return arma::dot(h, W*h)/(Ncum(M)*Ncum(M)*R*R*T*T*S*S);
}

//[[Rcpp::export]]
List fmvzeta3nc(const double& alpha,
                const arma::vec& beta,
                const int& R,
                const int& S,
                const int& T, 
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::mat& X1, 
                const arma::mat& X2, 
                const arma::mat& W,
                const int& Kx1, 
                const int& ninstr,
                const int& M, 
                const arma::vec& N, 
                const int& Pm,
                const arma::vec& Ncum){
  arma::vec X1b = X1*beta.head(Kx1);
  arma::mat dddG, ddG, dA, dZ, dddGX2;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Rabm(ninstr, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec X1bm = X1b.subvec(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm);
        Daym        += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          Rabm      += (dZ.t()*dA*arma::solve(Im - alpha*ddG, X1bm));
        }
      }
    }
    matM.col(m)    = Daym/(R*S) - Rabm/(R*S*T);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH3nc(const double& alpha,
                 const arma::vec& beta,
                 const int& R,
                 const int& S,
                 const int& T, 
                 List& distr, 
                 List& Ilist, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat dG, ddG, dddG, dA, dZ, dddGX2, ddA, ddAV;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm);
      for(int s(0); s < S; ++s){
        dG           = fGm(dnm, Nm);
        dA           = Im - alpha*dG;
        Daym        += (dZ.t()*dA*ym);
        Dgym        += (dZ.t()*dG*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddA        = Im - alpha*ddG;
          ddAV       = arma::solve(ddA, X1m);
          tmp1      += (dZ.t()*dG*ddAV);
          tmp2      += (dZ.t()*dA*arma::solve(ddA, ddG*ddAV));
          Ram       += (dZ.t()*dA*ddAV);
        } 
      }
    }
    matM.col(m)    = Daym/(R*S) - (Ram*beta)/(R*S*T);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + 1);
  deM.col(0)       = (T*Dgy - (tmp1 - tmp2)*beta)/(-Ncum(M)*R*S*T);
  deM.cols(1, Kx1) = Ra/(-Ncum(M)*R*S*T);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}
//*** Contextual, fixed effects
//[[Rcpp::export]]
arma::vec fbeta3fe(const double& alpha,
                   arma::vec& Day,
                   arma::mat& Ra,
                   const int& R,
                   const int& S,
                   const int& T, 
                   List& distr, 
                   List& Ilist, 
                   const arma::vec& y, 
                   const arma::mat& X1, 
                   const arma::mat& X2, 
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& Kx2, 
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ, dddGX2;
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm); dZ.each_row() -= mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm); dA.each_row() -= mean(dA, 0);
        Day         += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          Ra        += (dZ.t()*dA*arma::solve(Im - alpha*ddG, arma::join_rows(X1m, ddG*X2m)));
        }
      }
    }
  }
  Day     *= T;
  return arma::solve(Ra.t()*W*Ra, Ra.t()*(W*Day));
}

//[[Rcpp::export]]
double fgmm3fe(const double& alpha,
               const int& R,
               const int& S,
               const int& T, 
               List& distr, 
               List& Ilist, 
               const arma::vec& y, 
               const arma::mat& X1, 
               const arma::mat& X2, 
               const arma::mat& W,
               const int& Kx1, 
               const int& Kx2, 
               const int& ninstr,
               const int& M, 
               const arma::vec& N, 
               const int& Pm,
               const arma::vec& Ncum){
  arma::vec Day(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::vec beta = fbeta3fe(alpha, Day, Ra, R, S, T, distr, Ilist, y, X1, X2, W, Kx1, Kx2, M, N, Pm, Ncum);
  arma::vec h = Day - Ra*beta;
  return arma::dot(h, W*h)/(Ncum(M)*Ncum(M)*R*R*T*T*S*S);
}

//[[Rcpp::export]]
List fmvzeta3fe(const double& alpha,
                const arma::vec& beta,
                const int& R,
                const int& S,
                const int& T, 
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::mat& X1, 
                const arma::mat& X2, 
                const arma::mat& W,
                const int& Kx1, 
                const int& Kx2, 
                const int& ninstr,
                const int& M, 
                const arma::vec& N, 
                const int& Pm,
                const arma::vec& Ncum){
  arma::vec X1b = X1*beta.head(Kx1);
  arma::vec X2b = X2*beta.tail(Kx2);
  arma::mat dddG, ddG, dA, dZ, dddGX2;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Rabm(ninstr, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec X1bm = X1b.subvec(n1, n2);
    arma::vec X2bm = X2b.subvec(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm); dZ.each_row() -= mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm); dA.each_row() -= mean(dA, 0);
        Daym        += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          Rabm      += (dZ.t()*dA*arma::solve(Im - alpha*ddG, X1bm + ddG*X2bm));
        }
      }
    }
    matM.col(m)    = Daym/(R*S) - Rabm/(R*S*T);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH3fe(const double& alpha,
                 const arma::vec& beta,
                 const int& R,
                 const int& S,
                 const int& T, 
                 List& distr, 
                 List& Ilist, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& Kx2, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::mat dG, ddG, dddG, dA, dZ, dddGX2, ddA, ddAV, dGc;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1 + Kx2, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1 + Kx2, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1 + Kx2, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm); dZ.each_row() -= mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dG           = fGm(dnm, Nm);
        dGc          = dG.each_row() - mean(dG, 0);
        dA           = Im - alpha*dG; dA.each_row() -= mean(dA, 0);
        Daym        += (dZ.t()*dA*ym);
        Dgym        += (dZ.t()*dGc*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddA        = Im - alpha*ddG;
          ddAV       = arma::solve(ddA, arma::join_rows(X1m, ddG*X2m));
          tmp1      += (dZ.t()*dGc*ddAV);
          tmp2      += (dZ.t()*dA*arma::solve(ddA, ddG*ddAV));
          Ram       += (dZ.t()*dA*ddAV);
        }
      }
    }
    matM.col(m)    = Daym/(R*S) - (Ram*beta)/(R*S*T);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + Kx2 + 1);
  deM.col(0)             = (T*Dgy - (tmp1 - tmp2)*beta)/(-Ncum(M)*R*S*T);
  deM.cols(1, Kx1 + Kx2) = Ra/(-Ncum(M)*R*S*T);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}

//*** No Contextual, fixed effects
//[[Rcpp::export]]
arma::vec fbeta3ncfe(const double& alpha,
                     arma::vec& Day,
                     arma::mat& Ra,
                     const int& R,
                     const int& S,
                     const int& T, 
                     List& distr, 
                     List& Ilist, 
                     const arma::vec& y, 
                     const arma::mat& X1, 
                     const arma::mat& X2, 
                     const arma::mat& W,
                     const int& Kx1, 
                     const int& M, 
                     const arma::vec& N, 
                     const int& Pm,
                     const arma::vec& Ncum){
  arma::mat dddG, ddG, dA, dZ, dddGX2;
  for(int m(0); m < M; ++m){
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm); dZ.each_row() -= mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm); dA.each_row() -= mean(dA, 0);
        Day         += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          Ra        += (dZ.t()*dA*arma::solve(Im - alpha*ddG, X1m));
        }
      }
    }
  }
  Day     *= T;
  return arma::solve(Ra.t()*W*Ra, Ra.t()*(W*Day));
}

//[[Rcpp::export]]
double fgmm3ncfe(const double& alpha,
                 const int& R,
                 const int& S,
                 const int& T, 
                 List& distr, 
                 List& Ilist, 
                 const arma::vec& y, 
                 const arma::mat& X1, 
                 const arma::mat& X2, 
                 const arma::mat& W,
                 const int& Kx1, 
                 const int& ninstr,
                 const int& M, 
                 const arma::vec& N, 
                 const int& Pm,
                 const arma::vec& Ncum){
  arma::vec Day(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1, arma::fill::zeros);
  arma::vec beta = fbeta3ncfe(alpha, Day, Ra, R, S, T, distr, Ilist, y, X1, X2, W, Kx1, M, N, Pm, Ncum);
  arma::vec h = Day - Ra*beta;
  return arma::dot(h, W*h)/(Ncum(M)*Ncum(M)*R*R*T*T*S*S);
}

//[[Rcpp::export]]
List fmvzeta3ncfe(const double& alpha,
                  const arma::vec& beta,
                  const int& R,
                  const int& S,
                  const int& T, 
                  List& distr, 
                  List& Ilist, 
                  const arma::vec& y, 
                  const arma::mat& X1, 
                  const arma::mat& X2, 
                  const arma::mat& W,
                  const int& Kx1, 
                  const int& ninstr,
                  const int& M, 
                  const arma::vec& N, 
                  const int& Pm,
                  const arma::vec& Ncum){
  arma::vec X1b = X1*beta.head(Kx1);
  arma::mat dddG, ddG, dA, dZ, dddGX2;
  arma::mat matM(ninstr, M);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Rabm(ninstr, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    arma::vec X1bm = X1b.subvec(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2);  fZ(dZ, dddGX2, dddG, Pm); dZ.each_row() -= mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dA           = Im - alpha*fGm(dnm, Nm); dA.each_row() -= mean(dA, 0);
        Daym        += (dZ.t()*dA*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          Rabm      += (dZ.t()*dA*arma::solve(Im - alpha*ddG, X1bm));
        }
      }
    }
    matM.col(m)    = Daym/(R*S) - Rabm/(R*S*T);
  }
  return List::create(Named("sumM") = arma::sum(matM, 1), Named("sumMM") = matM*matM.t());
}

//[[Rcpp::export]]
List fmvzetaH3ncfe(const double& alpha,
                   const arma::vec& beta,
                   const int& R,
                   const int& S,
                   const int& T, 
                   List& distr, 
                   List& Ilist, 
                   const arma::vec& y, 
                   const arma::mat& X1, 
                   const arma::mat& X2, 
                   const arma::mat& W,
                   const int& Kx1, 
                   const int& ninstr,
                   const int& M, 
                   const arma::vec& N, 
                   const int& Pm,
                   const arma::vec& Ncum){
  arma::mat dG, ddG, dddG, dA, dZ, dddGX2, ddA, ddAV, dGc;
  arma::mat matM(ninstr, M);
  arma::vec Dgy(ninstr, arma::fill::zeros);
  arma::mat Ra(ninstr, Kx1, arma::fill::zeros);
  arma::mat tmp1(ninstr, Kx1, arma::fill::zeros);
  arma::mat tmp2(ninstr, Kx1, arma::fill::zeros);
  for(int m(0); m < M; ++m){
    arma::vec Daym(ninstr, arma::fill::zeros);
    arma::vec Dgym(ninstr, arma::fill::zeros);
    arma::mat Ram(ninstr, Kx1, arma::fill::zeros);
    int Nm         = N(m);
    int n1         = Ncum(m);
    int n2         = Ncum(m + 1) - 1;
    arma::mat dnm  = distr(m);
    arma::mat Im   = Ilist(m);
    arma::vec ym   = y.subvec(n1, n2);
    arma::mat X1m  = X1.rows(n1, n2);
    arma::mat X2m  = X2.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dddG           = fGm(dnm, Nm);
      dddGX2         = dddG*X2m;
      dZ             = arma::join_rows(X1m, dddGX2); fZ(dZ, dddGX2, dddG, Pm); dZ.each_row() -= mean(dZ, 0);
      for(int s(0); s < S; ++s){
        dG           = fGm(dnm, Nm);
        dGc          = dG.each_row() - mean(dG, 0);
        dA           = Im - alpha*dG; dA.each_row() -= mean(dA, 0);
        Daym        += (dZ.t()*dA*ym);
        Dgym        += (dZ.t()*dGc*ym);
        for(int t(0); t < T; ++t){
          ddG        = fGm(dnm, Nm);
          ddA        = Im - alpha*ddG;
          ddAV       = arma::solve(ddA, X1m);
          tmp1      += (dZ.t()*dGc*ddAV);
          tmp2      += (dZ.t()*dA*arma::solve(ddA, ddG*ddAV));
          Ram       += (dZ.t()*dA*ddAV);
        }
      }
    }
    matM.col(m)    = Daym/(R*S) - (Ram*beta)/(R*S*T);
    Dgy           += Dgym;
    Ra            += Ram;
  }
  arma::mat deM(ninstr, Kx1 + 1);
  deM.col(0)       = (T*Dgy - (tmp1 - tmp2)*beta)/(-Ncum(M)*R*S*T);
  deM.cols(1, Kx1) = Ra/(-Ncum(M)*R*S*T);
  return List::create(Named("sumM")  = arma::sum(matM, 1), 
                      Named("sumMM") = matM*matM.t(),
                      Named("derM")  = deM);
}