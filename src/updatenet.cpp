// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#define NDEBUG 1
#include <RcppEigen.h>
#include "tools.h"

using namespace Rcpp;
//using namespace arma;
using namespace std;
//using namespace Eigen;
//using namespace Numer;

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
    Eigen::MatrixXd eyeM = Eigen::MatrixXd::Identity(Nm,Nm);
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
      
      arma::vec IndexMJ = IndexMIJ[ci];
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
        arma::vec index_col = IndexMJ.subvec(Blocks(cj) + 1, Blocks(cj+1));
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
      
      arma::vec IndexMJ = IndexMIJ[ci];
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
        arma::vec index_col = IndexMJ.subvec(Blocks(cj) + 1, Blocks(cj+1));
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