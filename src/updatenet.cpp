// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#define NDEBUG 1
#include <RcppEigen.h>
#include "tools.h"
#include "vMF.h"

using namespace Rcpp;
//using namespace arma;
using namespace std;
//using namespace Eigen;
//using namespace Numer;
//In this file there are the functions which update the network and the prior distribution of the network


// compute the block of the line
arma::vec cBlock(
    const int& nupmax, 
    const int& NJ, 
    int& NJeff) {
  arma::vec out = cumsum(Rcpp::RcppArmadillo::sample(arma::linspace(1,nupmax,nupmax), NJ, true)) - 1;
  out = out.elem(arma::find(out <= (NJ - 1)));
  out = arma::join_cols(-arma::ones(1), out);
  NJeff = out.n_elem - 1;
  if(out(NJeff) != (NJ - 1)){
    out = arma::join_cols(out, arma::ones(1)*(NJ - 1));
    NJeff += 1;
  }
  
  return out;
}

// USE A WAY TO QUICLY COMPUTE THE DET AND THE INVERSE 
// A Structural Model for the Coevolution of Networks and Behavior

// Return one G from G|y one step of MCMC
// With contextual
void updGnorm (List& Gnorm,
               const List& prior,
               const List& ListIndex,
               const Rcpp::IntegerVector& N,
               const int& M,
               const List& y,
               List& A,
               List& Ay,
               const List& Xb,
               const List& Xgamma,
               const double& alpha,
               const double& sigma2) {
  for(int m(0); m<M; ++m){
    int Nm = N(m);
    Eigen::MatrixXd eyeM = Eigen::MatrixXd::Identity(Nm,Nm);
    Eigen::VectorXd ym = y[m];
    Eigen::VectorXd Xmb = Xb[m];  
    Eigen::VectorXd Xmgamma = Xgamma[m];
    Eigen::MatrixXd priorm = prior(m);
    
    //Copy G[m]
    
    Eigen::MatrixXd Gmnorm = Gnorm(m);
    List IndexM = ListIndex[m];
    arma::vec IndexMI = IndexM[0]; // index row to update
    List IndexMIJ = IndexM[1]; // in index tow i, the colums to update
    int NI = IndexMI.n_elem;
    if (NI == 0) continue;
    
    Eigen::MatrixXd Am       = eyeM-alpha*Gmnorm;
    Eigen::MatrixXd invAm    = Am.inverse();
    double detAm             = Am.determinant();
    
    for(int ci(0); ci<NI; ++ci){
      int i = IndexMI(ci);
      Eigen::RowVectorXd Gmrowi = Gmnorm.row(i).array().ceil();
      int ni                    = Gmrowi.sum();
      arma::vec IndexMJ         = IndexMIJ[ci];
      int NJ                    = IndexMJ.n_elem;
      
      for(int cj(0); cj<NJ; ++cj){
        int j    = IndexMJ(cj);
        double l = Am(i, j);
        int gij1 = (l == 0); // new link
        int ni1  =  ni + gij1 - (gij1 == 0);
        double Aii1 = 1, Aij1 = -alpha*gij1*1.0;
        if (ni != 0 && ni1 != 0) {
          Aii1  = ni1*1.0/ni;
          Aij1  = -alpha*gij1*1.0/ni;
        }
        double eps_1 = Aii1 - 1.0;
        double eps_2 = Aij1 - l;
        
        Eigen::MatrixXd Am1       = Am;
        Am1(i, i)                 = Aii1;
        Am1(i, j)                 = Aij1;
        Eigen::MatrixXd invAm1    = invmodijk (invAm, i, i, j, eps_1, eps_2);
        double detAm1             = detmodijk (detAm, invAm, i, i, j, eps_1, eps_2);
        
        // normalization
        Am1.row(i)               /= Aii1;
        invAm1.col(i)            *= Aii1; 
        detAm1                   /= Aii1;
        
        Eigen::MatrixXd Gmnorm1   = Gmnorm;
        Gmnorm1(i,j)              = gij1;
        Gmnorm1.row(i)            = Gmnorm1.row(i).array().ceil();
        double rs                 = Gmnorm1.row(i).sum();
        if(rs>0){
          Gmnorm1.row(i)         /= rs;
        }
        
        //cout<<ci<<"/"<<NI<<" "<<cj<<"/"<<NJ<<endl;
        Eigen::Array2d PROB;
        Eigen::VectorXd tmpxb = Am*ym - Xmb - Gmnorm*Xmgamma;
        double tmpdot = (double) (0.5/sigma2)*tmpxb.transpose()*tmpxb; 
        PROB(1 - gij1)        = -0.5*Nm*log(sigma2) + log(detAm) - tmpdot;
        
        
        tmpxb                 = Am1*ym - Xmb - Gmnorm1*Xmgamma;
        tmpdot                =  (double) (0.5/sigma2)*tmpxb.transpose()*tmpxb; 
        PROB(gij1)            = -0.5*Nm*log(sigma2) + log(detAm1) - tmpdot;
        
        //Normalization and prior
        double maxprob = PROB.maxCoeff(); 
        PROB -= maxprob;
        PROB = PROB.exp();
        PROB(0) *= (1-priorm(i,j));
        PROB(1) *= priorm(i,j);
        double sumprob = PROB.sum();
        PROB /= sumprob;
        
        int gijsel = (unif_rand()<PROB(1));//arma::conv_to<int>::from(Rcpp::RcppArmadillo::sample_main(arma::linspace(0,1,2),1,false,PROB));
        if (gijsel == gij1) {
          Am     = Am1;
          invAm  = invAm1; 
          detAm  = detAm1;
          Gmnorm = Gmnorm1;
          ni     = ni1;
        }
      }
    }
    //Gnorm
    Gnorm[m] = Gmnorm;
    
    // A
    A[m] = Am;
    
    // Ay
    Eigen::VectorXd Amy = Am*ym;
    Ay[m] = Amy;
    
  }
}
// Without contextual
void updGnormnoc (List& Gnorm,
                  const List& prior,
                  const List& ListIndex,
                  const Rcpp::IntegerVector& N,
                  const int& M,
                  const List& y,
                  List& A,
                  List& Ay,
                  const List& Xb,
                  const double& alpha,
                  const double& sigma2) {
  for(int m(0); m<M; ++m){
    int Nm = N(m);
    Eigen::MatrixXd eyeM = Eigen::MatrixXd::Identity(Nm,Nm);
    Eigen::VectorXd ym = y[m];
    Eigen::VectorXd Xmb = Xb[m];  
    Eigen::MatrixXd priorm = prior(m);
    
    //Copy G[m]
    
    Eigen::MatrixXd Gmnorm = Gnorm(m);
    List IndexM = ListIndex[m];
    arma::vec IndexMI = IndexM[0]; // index row to update
    List IndexMIJ = IndexM[1]; // in index tow i, the colums to update
    int NI = IndexMI.n_elem;
    if (NI == 0) continue;
    
    Eigen::MatrixXd Am       = eyeM-alpha*Gmnorm;
    Eigen::MatrixXd invAm    = Am.inverse();
    double detAm             = Am.determinant();
    
    for(int ci(0); ci<NI; ++ci){
      int i = IndexMI(ci);
      Eigen::RowVectorXd Gmrowi = Gmnorm.row(i).array().ceil();
      int ni                    = Gmrowi.sum();
      arma::vec IndexMJ         = IndexMIJ[ci];
      int NJ                    = IndexMJ.n_elem;
      
      for(int cj(0); cj<NJ; ++cj){
        int j    = IndexMJ(cj);
        double l = Am(i, j);
        int gij1 = (l == 0); // new link
        int ni1  =  ni + gij1 - (gij1 == 0);
        double Aii1 = 1, Aij1 = -alpha*gij1*1.0;
        if (ni != 0 && ni1 != 0) {
          Aii1  = ni1*1.0/ni;
          Aij1  = -alpha*gij1*1.0/ni;
        }
        double eps_1 = Aii1 - 1.0;
        double eps_2 = Aij1 - l;
        
        Eigen::MatrixXd Am1       = Am;
        Am1(i, i)                 = Aii1;
        Am1(i, j)                 = Aij1;
        Eigen::MatrixXd invAm1    = invmodijk (invAm, i, i, j, eps_1, eps_2);
        double detAm1             = detmodijk (detAm, invAm, i, i, j, eps_1, eps_2);
        
        // normalization
        Am1.row(i)               /= Aii1;
        invAm1.col(i)            *= Aii1; 
        detAm1                   /= Aii1;
        
        Eigen::MatrixXd Gmnorm1   = Gmnorm;
        Gmnorm1(i,j)              = gij1;
        Gmnorm1.row(i)            = Gmnorm1.row(i).array().ceil();
        double rs                 = Gmnorm1.row(i).sum();
        if(rs>0){
          Gmnorm1.row(i)         /= rs;
        }
        
        //cout<<ci<<"/"<<NI<<" "<<cj<<"/"<<NJ<<endl;
        Eigen::Array2d PROB;
        Eigen::VectorXd tmpxb = Am*ym - Xmb;
        double tmpdot = (double) (0.5/sigma2)*tmpxb.transpose()*tmpxb; 
        PROB(1 - gij1)        = -0.5*Nm*log(sigma2) + log(detAm) - tmpdot;
        
        
        tmpxb                 = Am1*ym - Xmb;
        tmpdot                =  (double) (0.5/sigma2)*tmpxb.transpose()*tmpxb; 
        PROB(gij1)            = -0.5*Nm*log(sigma2) + log(detAm1) - tmpdot;
        
        //Normalization and prior
        double maxprob = PROB.maxCoeff(); 
        PROB -= maxprob;
        PROB = PROB.exp();
        PROB(0) *= (1-priorm(i,j));
        PROB(1) *= priorm(i,j);
        double sumprob = PROB.sum();
        PROB /= sumprob;
        
        int gijsel = (unif_rand()<PROB(1));//arma::conv_to<int>::from(Rcpp::RcppArmadillo::sample_main(arma::linspace(0,1,2),1,false,PROB));
        if (gijsel == gij1) {
          Am     = Am1;
          invAm  = invAm1; 
          detAm  = detAm1;
          Gmnorm = Gmnorm1;
          ni     = ni1;
        }
      }
    }
    //Gnorm
    Gnorm[m] = Gmnorm;
    
    // A
    A[m] = Am;
    
    // Ay
    Eigen::VectorXd Amy = Am*ym;
    Ay[m] = Amy;
    
  }
}
// In block  withou contextual
void updGnormblock (List& Gnorm,
                    const List& prior,
                    const List& ListIndex,
                    const Rcpp::IntegerVector& N,
                    const int& M,
                    const List& y,
                    List& A,
                    List& Ay,
                    const List& Xb,
                    const List& Xgamma,
                    const double& alpha,
                    const double& sigma2,
                    const int& nupmax) {
  for(int m(0); m<M; ++m){
    int Nm = N(m);
    Eigen::MatrixXd eyeM = Eigen::MatrixXd::Identity(Nm, Nm);
    Eigen::VectorXd ym = y[m];
    Eigen::VectorXd Xmb = Xb[m];  
    Eigen::VectorXd Xmgamma = Xgamma[m];
    Eigen::MatrixXd priorm = prior(m);
    
    //Copy G[m]
    
    Eigen::MatrixXd Gmnorm = Gnorm(m);
    List IndexM = ListIndex[m];
    arma::vec IndexMI = IndexM[0];
    List IndexMIJ = IndexM[1];
    
    int NI = IndexMI.n_elem;
    for(int ci(0); ci<NI; ++ci){ // line
      int i = IndexMI(ci);
      // compute N cofactor for row i
      Eigen::VectorXd Cofator(N(m));
      Eigen::MatrixXd Am = eyeM-alpha*Gmnorm;
      removeRow(Am,i);
      
      arma::Col<int> IndexMJ = IndexMIJ[ci];
      int NJ = IndexMJ.n_elem;
      
      for(int j(0); j < N(m); ++j) {
        Eigen::MatrixXd Amtemp = Am;
        removeColumn(Amtemp,j);
        Cofator(j) = Amtemp.determinant()*pow(-1,i+j);
      }
      
      // permutations and blocks
      int NJeff;
      IndexMJ = Rcpp::RcppArmadillo::sample(IndexMJ, NJ, false);
      arma::vec Blocks  = cBlock(nupmax,  NJ,  NJeff);
      
      
      for(int cj(0); cj<NJeff; ++cj){
        arma::Col<int> index_col = IndexMJ.subvec(Blocks(cj) + 1, Blocks(cj+1));
        int nupdate         = Blocks(cj+1) - Blocks(cj);
        int poss            = pow(2, nupdate);          // number of possibility
        arma::mat cposs     = possentries(nupdate, poss); // the possibility 
        
        arma::vec PROB(poss), prior_block(poss);
        
        for (int l(0); l < poss; ++ l) {
          double prior_blockl = 1;
          updselel(Gmnorm, prior_blockl, priorm, i, index_col, cposs.row(l), nupdate);
          
          Gmnorm.row(i)=Gmnorm.row(i).array().ceil();
          double rs = Gmnorm.row(i).sum();
          if(rs>0){
            Gmnorm.row(i) /= rs;
          }
          
          Eigen::MatrixXd Am = eyeM-alpha*Gmnorm;
          Eigen::VectorXd Amy = Am*ym;
          Eigen::VectorXd GmXmgamma = Gmnorm*Xmgamma;
          PROB(l) =  ((Am.row(i)*Cofator).log()  - (0.5/sigma2)*(Amy.transpose()*(Amy - 2*Xmb - 2*GmXmgamma) + GmXmgamma.transpose()*(2*Xmb + GmXmgamma)))(0);
          prior_block(l) = prior_blockl;
        }
        
        //Normalization and prior
        PROB = exp(PROB - max(PROB));
        PROB = PROB%prior_block;
        PROB = PROB / sum(PROB);
        
        
        int l = sum(Rcpp::RcppArmadillo::sample_main(arma::linspace(0, poss - 1, poss), 1, true, PROB));
        double prior_blockl = 1;
        updselel(Gmnorm, prior_blockl, priorm, i, index_col, cposs.row(l), nupdate);
        
        Gmnorm.row(i)=Gmnorm.row(i).array().ceil();
        double rs = Gmnorm.row(i).sum();
        if(rs>0){
          Gmnorm.row(i) /= rs;
        }       
      }
    }
    //Gnorm
    Gnorm[m] = Gmnorm;
    
    // A
    Eigen::MatrixXd Am = eyeM-alpha*Gmnorm;
    A[m] = Am;
    
    // Ay
    Eigen::VectorXd Amy = Am*ym;
    Ay[m] = Amy;
  }
}
// in block without contextual
void updGnormblocknoc (List& Gnorm,
                       const List& prior,
                       const List& ListIndex,
                       const Rcpp::IntegerVector& N,
                       const int& M,
                       const List& y,
                       List& A,
                       List& Ay,
                       const List& Xb,
                       const double& alpha,
                       const double& sigma2,
                       const int& nupmax) {
  for(int m(0); m<M; ++m){
    int Nm = N(m);
    Eigen::MatrixXd eyeM = Eigen::MatrixXd::Identity(Nm,Nm);
    Eigen::VectorXd ym = y[m];
    Eigen::VectorXd Xmb = Xb[m];  
    Eigen::MatrixXd priorm = prior(m);
    
    //Copy G[m]
    
    Eigen::MatrixXd Gmnorm = Gnorm(m);
    List IndexM = ListIndex[m];
    arma::vec IndexMI = IndexM[0];
    List IndexMIJ = IndexM[1];
    
    int NI = IndexMI.n_elem;
    for(int ci(0); ci<NI; ++ci){ // line
      int i = IndexMI(ci);
      // compute N cofactor for row i
      Eigen::VectorXd Cofator(N(m));
      Eigen::MatrixXd Am = eyeM-alpha*Gmnorm;
      removeRow(Am,i);
      
      arma::Col<int> IndexMJ = IndexMIJ[ci];
      int NJ = IndexMJ.n_elem;
      
      for(int j(0); j < N(m); ++j) {
        Eigen::MatrixXd Amtemp = Am;
        removeColumn(Amtemp,j);
        Cofator(j) = Amtemp.determinant()*pow(-1,i+j);
      }
      
      // permutations and blocks
      int NJeff;
      IndexMJ = Rcpp::RcppArmadillo::sample(IndexMJ, NJ, false);
      arma::vec Blocks  = cBlock(nupmax,  NJ,  NJeff);
      
      
      for(int cj(0); cj<NJeff; ++cj){
        arma::Col<int> index_col = IndexMJ.subvec(Blocks(cj) + 1, Blocks(cj+1));
        int nupdate         = Blocks(cj+1) - Blocks(cj);
        int poss            = pow(2, nupdate);          // number of possibility
        arma::mat cposs     = possentries(nupdate, poss); // the possibility 
        
        arma::vec PROB(poss), prior_block(poss);
        
        for (int l(0); l < poss; ++ l) {
          double prior_blockl = 1;
          updselel(Gmnorm, prior_blockl, priorm, i, index_col, cposs.row(l), nupdate);
          
          Gmnorm.row(i)=Gmnorm.row(i).array().ceil();
          double rs = Gmnorm.row(i).sum();
          if(rs>0){
            Gmnorm.row(i) /= rs;
          }
          
          Eigen::MatrixXd Am = eyeM-alpha*Gmnorm;
          Eigen::VectorXd Amy = Am*ym;
          PROB(l) =  ((Am.row(i)*Cofator).log()  - (0.5/sigma2)*(Amy.transpose()*(Amy - 2*Xmb)))(0);
          prior_block(l) = prior_blockl;
        }
        
        //Normalization and prior
        PROB = exp(PROB - max(PROB));
        PROB = PROB%prior_block;
        PROB = PROB / sum(PROB);
        
        
        int l = sum(Rcpp::RcppArmadillo::sample_main(arma::linspace(0, poss - 1, poss), 1, true, PROB));
        double prior_blockl = 1;
        updselel(Gmnorm, prior_blockl, priorm, i, index_col, cposs.row(l), nupdate);
        
        Gmnorm.row(i)=Gmnorm.row(i).array().ceil();
        double rs = Gmnorm.row(i).sum();
        if(rs>0){
          Gmnorm.row(i) /= rs;
        }       
      }
    }
    //Gnorm
    Gnorm[m] = Gmnorm;
    
    // A
    Eigen::MatrixXd Am = eyeM-alpha*Gmnorm;
    A[m] = Am;
    
    // Ay
    Eigen::VectorXd Amy = Am*ym;
    Ay[m] = Amy;
  }
}


// Update rho using probit or logit
void updrhopl(List& Gnorm,
              List& prior,
              List& G0obs,
              List& ListIndex,
              arma::vec& rho,
              Eigen::VectorXd& lFdZrhoE1,
              Eigen::VectorXd& lFdZrhoE0,
              const Rcpp::NumericVector weight,
              const arma::mat& dZ,
              const arma::vec& murho,
              const arma::mat& iVrho,
              const arma::mat& jumprho,
              const int& Krho,
              const Rcpp::IntegerVector& N,
              const int& M,
              double& rhoaccept,
              const int& type,
              const bool& Afixed,
              const Eigen::ArrayXd& G0obsvec) {
  Eigen::VectorXd netA           = frMtoV(Gnorm,  N, M);
  arma::vec rhost                = Fmvnorm(Krho, rho, jumprho);
  arma::vec dZrhost              = dZ*rhost;
  NumericVector dZrhostR         = wrap(dZrhost);
  NumericVector lFdZrhostR1, lFdZrhostR0;
  
  if (type == 1) { //probit
    lFdZrhostR1    = Rcpp::pnorm(dZrhostR, 0, 1, true, true);
    lFdZrhostR0    = Rcpp::pnorm(dZrhostR, 0, 1, false, true);
  } else { // is supposed to be 2: logit
    lFdZrhostR1    = Rcpp::plogis(dZrhostR, 0, 1, true, true);
    lFdZrhostR0    = Rcpp::plogis(dZrhostR, 0, 1, false, true);
  }
  lFdZrhostR1                    = lFdZrhostR1*weight;
  lFdZrhostR0                    = lFdZrhostR0*weight;
  Eigen::VectorXd lFdZrhostE1    = as<Eigen::VectorXd>(lFdZrhostR1);
  Eigen::VectorXd lFdZrhostE0    = as<Eigen::VectorXd>(lFdZrhostR0);
  // cout<<weight(1)<<endl;
  // acceptance rate
  Eigen::VectorXd tmp = (netA.array() > 0).select(lFdZrhostE1 - lFdZrhoE1, lFdZrhostE0 - lFdZrhoE0);
  if(Afixed){
    tmp               = (G0obsvec == 0).select(tmp, 0);
  }
  double lalpharho1   = tmp.sum();
  double lalpharho2   = 0.5*(arma::dot(rho - murho, iVrho*(rho - murho)) - arma::dot(rhost - murho, iVrho*(rhost - murho)));
  double logalpharho  = lalpharho1 + lalpharho2;
  
  if(unif_rand()<exp(logalpharho)){
    prior             = frVtoM(lFdZrhostE1.array().exp(),  N,  M);
    rho               = rhost;
    lFdZrhoE1         = lFdZrhostE1;
    lFdZrhoE0         = lFdZrhostE0;
    ListIndex         = fListIndex(prior, G0obs, M,  N);
    rhoaccept        += 1;     //Increase acceptance number to 1 
  }
  // cout<<rho.t()<<endl;
  // arma::mat prior1 = prior[1];
  // cout<<prior1.row(1)<<endl;
}


// Update rho using the latent space model
// ARD are observed in the whole population
void updrhoARD(List& Gnorm,
               List& prior,
               List& G0obs,
               List& ListIndex,
               List& rho,
               const List& d,
               const arma::vec& zeta,
               const List& murho,
               const List& iVrho,
               const List& jumprho,
               const arma::vec& Krho,
               List& neighbor,
               List& weight,
               List& iARD,
               List& inonARD,
               const Rcpp::IntegerVector& N,
               const Rcpp::IntegerVector& N1,
               const int& M,
               const Rcpp::IntegerVector& P,
               arma::vec& rhoaccept,
               const arma::vec& type) {
  
  for(int m(0); m < M; ++m) {
    int N1m              = N1(m);
    int Nm               = N(m);
    int N2m              = Nm - N1m;
    arma::vec rhom       = rho[m];
    arma::mat jumprhom   = jumprho[m];
    arma::vec rhomst     = Fmvnorm(Krho(m), rhom, jumprhom);
    
    arma::mat priorm     = prior[m];
    arma::mat Gm         = Gnorm[m];
    Gm                   = arma::ceil(Gm);
    arma::mat G0obsm     = G0obs[m];
    arma::vec murhom     = murho[m];
    arma::mat iVrhom     = iVrho[m];
    
    arma::umat neighborm = neighbor[m]; 
    arma::mat weightm    = weight[m];
    arma::uvec iARDm     = iARD[m];
    arma::uvec inonARDm  = inonARD[m];
    
    arma::mat priorstm;
    if(type(m) == 0) { //d and zeta vary
      double zetam      = exp(rhomst(0));
      double logCpzetam = logCpvMFcpp(P(m), zetam);
      
      arma::vec num     = rhomst.tail(N1m);
      arma::vec dm      = (Nm*1.0/N1m)*exp(num - logCpzetam)*sum(exp(num));
      
      arma::mat zm      = rhomst.subvec(1, N1m*P(m));
      zm.reshape(N1m, P(m));

      priorstm          =  fdnetARD(zm, num, dm, N1m, N2m, Nm, P(m), zetam, logCpzetam, neighborm, weightm, iARDm, inonARDm);
    }
    if(type(m) == 1) { //d fixed and zeta varies
      double zetam      = exp(rhomst(0));
      double logCpzetam = logCpvMFcpp(P(m), zetam);
      
      arma::vec dm      = d[m];
      arma::vec num     = log(dm) + 0.5*logCpzetam + 0.5*log(N1m*1.0/Nm) - 0.5*log(sum(dm));
      
      arma::mat zm      = rhomst.subvec(1, N1m*P(m));
      zm.reshape(N1m, P(m));
      
      priorstm          =  fdnetARD(zm, num, dm, N1m, N2m, Nm, P(m), zetam, logCpzetam, neighborm, weightm, iARDm, inonARDm);
    }
    if(type(m) == 2) { //d varies and zeta fixed
      double zetam      = zeta(m);
      double logCpzetam = logCpvMFcpp(P(m), zetam);
      
      arma::vec num     = rhomst.tail(N1m);
      arma::vec dm      = (Nm*1.0/N1m)*exp(num - logCpzetam)*sum(exp(num));
      
      arma::mat zm      = rhomst.subvec(0, N1m*P(m) - 1);
      zm.reshape(N1m, P(m));
      
      priorstm          =  fdnetARD(zm, num, dm, N1m, N2m, Nm, P(m), zetam, logCpzetam, neighborm, weightm, iARDm, inonARDm);
    }
    if(type(m) == 3) { //d and zeta fixed
      double zetam      = zeta(m);
      double logCpzetam = logCpvMFcpp(P(m), zetam);
      
      arma::vec dm      = d[m];
      arma::vec num     = log(dm) + 0.5*logCpzetam + 0.5*log(N1m*1.0/Nm) - 0.5*log(sum(dm));
      
      arma::mat zm      = rhomst;
      zm.reshape(N1m, P(m));
      
      priorstm          =  fdnetARD(zm, num, dm, N1m, N2m, Nm, P(m), zetam, logCpzetam, neighborm, weightm, iARDm, inonARDm);
    }
    
    double lalpharho1   = arma::accu((log(priorstm%(2*Gm - 1) + 1 - Gm) - log(priorm%(2*Gm - 1) + 1 - Gm))%(1 - G0obsm));
    double lalpharho2   = 0.5*(arma::dot(rhom - murhom, iVrhom*(rhom - murhom)) - arma::dot(rhomst - murhom, iVrhom*(rhomst - murhom)));
    
    double logalpharho  = lalpharho1 + lalpharho2;
    
    if(unif_rand()<exp(logalpharho)){
      prior[m]          = priorstm;
      rho[m]            = rhomst;
      rhoaccept(m)     += 1;     //Increase acceptance number to 1
    }
  }
  ListIndex             = fListIndex(prior, G0obs, M,  N);
}
