// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>

using namespace Rcpp;
//using namespace arma;
using namespace std;
//using namespace Eigen;
//using namespace Numer;

// USE A WAY TO QUICLY COMPUTE THE DET AND THE INVERSE 
// Hsieh et al. (2019), A Structural Model for the Coevolution of Networks and Behavior



//typedef Eigen::Map<Eigen::MatrixXd> MapMat;
//typedef Eigen::Map<Eigen::VectorXd> MapVec;

// gibbstep draws one graph from G|y X THETA distribution
// G is the matrix to update
// prior is the prior matrix with pij = prior probability that Gij = 1
// theta is the parameter vector
// y is the outcome
// X contains the explanatory variables with column 1 = intercept
// intercept is a boolean variable = 1 if X contains intercept
// c_index is the index of X that will be usen in contextual variable

//Fast multivariate normal sampling
arma::vec Fmvnorm(const double& dim, arma::vec u, arma::mat sigma) {
  arma::vec x = arma::randn(dim,1);
  return arma::chol(sigma).t()*x + u;
}

// Update invV if (i, j) is modified
Eigen::MatrixXd invmodij (const Eigen::MatrixXd& invM,
                          const unsigned int& i,
                          const unsigned int& j,
                          const double& eps_) {
  return invM - eps_*invM.col(i)*invM.row(j)/(1. + eps_*invM(j,i));
}

// Update invV if (i, j) and (i, k) are modified
Eigen::MatrixXd invmodijk (const Eigen::MatrixXd& invM,
                           const unsigned int& i,
                           const unsigned int& j,
                           const unsigned int& k,
                           const double& eps_1,
                           const double& eps_2) {
  Eigen::MatrixXd tmp = invM - eps_1*invM.col(i)*invM.row(j)/(1 + eps_1*invM(j,i));
  return tmp - eps_2*tmp.col(i)*tmp.row(k)/(1 + eps_2*tmp(k,i));
}

// Update detM if (i, j) is modified
double detmodij (const double& detM,
                 const Eigen::MatrixXd& invM,
                 const unsigned int& i,
                 const unsigned int& j,
                 const double& eps_) {
  return (1 + eps_*invM(j,i))*detM;
}

// Update detM if (i, j) and (i, k) are modified
double detmodijk (const double& detM,
                  const Eigen::MatrixXd& invM,
                  const unsigned int& i,
                  const unsigned int& j,
                  const unsigned int& k,
                  const double& eps_1,
                  const double& eps_2) {
  return (1 + eps_2*invM(k,i) - eps_1*eps_2*invM(k,i)*invM(j,i))*(1 + eps_1*invM(j,i))*detM;
}


// remove row and col from Eigen matrix
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();
  
  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;
  
  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}

// given a number of entries, 
// return all the possibilities of entries
arma::mat possentries(int& nupdate, const int& pow_nupdate) 
{ 
  /*nupdate of power set of a set with nupdate 
   n is (2**n -1)*/
  int counter, j; 
  
  arma::mat output(pow_nupdate, nupdate);
  
  
  for(counter = 0; counter < pow_nupdate; counter++) 
  { 
    arma::uvec temp(nupdate);
    int i = 0;
    
    for(j = 0; j < nupdate; j++) 
    { 
      /* Check if jth bit in the counter is set 
       If set then print jth element from set */
      if(counter & (1<<j)){
        temp(i) = j;
        i ++;
      } 
    }
    arma::rowvec outputc = arma::zeros<arma::rowvec>(nupdate); 
    outputc.elem(temp.head(i))  = arma::ones<arma::rowvec>(i);
    output.row(counter) = outputc;
  }
  
  return output;
} 

// update ligne i and col in index_col by newvalue
void updselel (
    Eigen::MatrixXd& Gmnorm,
    double& prior_blockl,
    const Eigen::MatrixXd priorm,
    const int& i,
    const arma::vec& index_col,
    const arma::rowvec& newval,
    const int& nupdate
){
  
  for (int j(0); j < nupdate; ++ j) {
    Gmnorm(i, index_col(j)) = newval(j);
    prior_blockl *= (priorm(i, index_col(j))*newval(j) + (1-priorm(i, index_col(j)))*(1-newval(j)));
  }
}

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
               const double& sigma2,
               const double& kbeta,
               const double& kv,
               const int t) {
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
                    const double& kbeta,
                    const double& kv,
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
          PROB(l) =  ((Am.row(i)*Cofator).log()  - (0.5/sigma2)*((Amy.transpose()*(Amy - 2*Xmb - 2*GmXmgamma) + GmXmgamma.transpose()*(2*Xmb + GmXmgamma))))(0);
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

// update theta
void updtheta (arma::vec& theta,
               List& Vtheta,
               List& Xb,
               List& Xgamma,
               const double& sigma2,
               const double& kv,
               const double& kbeta,
               const List& Xone,
               const List& X,
               const List& Ay,
               const List& V,
               const arma::vec& invsigmathetatheta0,
               const arma::mat& invsigmatheta,
               const double& M) {
  arma::vec Amy = Ay[0];
  arma::mat Vm  = V[0];
  arma::vec VAY = Vm.t()*Amy;
  arma::mat VV = Vm.t()*Vm;
  
  for (int m(1); m<M; ++m) {
    arma::vec Amy = Ay[m];
    arma::mat Vm  = V[m];
    VAY += Vm.t()*Amy;
    VV += Vm.t()*Vm;
  }
  
  const arma::vec muhtheta = VAY + invsigmathetatheta0;
  const arma::mat sigmahtheta = arma::inv(VV + invsigmatheta);
  theta = Fmvnorm(kv,sigmahtheta*muhtheta,sigma2*sigmahtheta);
  
  for(int m(0); m<M; ++m){
    arma::mat Vm = V[m];
    arma::mat Xonem = Xone[m];
    arma::mat Xm = X[m];
    arma::vec Vmtheta = Vm*theta;
    
    // Vtheta
    Vtheta[m] = Vmtheta;
    
    //Xb
    Xb[m] = Xonem*theta.head(kbeta);
    
    //Xgamma
    Xgamma[m] = Xm*theta.tail(kv-kbeta);
  }
}


// update sigma
void updsigma2 (double& sigma2,
                const arma::vec& theta,
                const double& a,
                const double& b,
                const arma::vec theta0,
                const arma::mat& invsigmatheta,
                const List& Ay,
                const List& Vtheta,
                const double& sumN,
                const double& M) {
  // ah
  double ah = 0.5*(a+sumN);
  arma::vec temp1 = theta - theta0;
  
  // bh
  double bh = b + arma::dot(temp1, invsigmatheta*temp1);
  
  arma::vec Amy = Ay[0];
  arma::vec Vmtheta  = Vtheta[0];
  arma::vec temp2 = Amy - Vmtheta;
  bh += arma::dot(temp2, temp2);
  for (int m(1); m<M; ++m) {
    arma::vec Amy = Ay[m];
    arma::vec Vmtheta  = Vtheta[m];
    temp2 = Amy - Vmtheta;
    bh += arma::dot(temp2, temp2);
  }
  bh *= 0.5;
  
  // draw sigma
  NumericVector sigma2temps = rgamma(1,ah,1/bh);
  sigma2 = 1/sigma2temps(0);
}


// update zeta
void updzeta (double& zeta,
              double& alpha,
              List& A,
              double& sumlogdetA,
              List& Ay,
              const List& Gnorm,
              const List& y,
              const double& sigma2,
              const List& Vtheta,
              const double& jumpzeta,
              double& zetaaccept,
              const double& zeta0,
              const double& invsigma2zeta,
              const Rcpp::IntegerVector N,
              const double M) {
  
  double zetastart =  R::rnorm(zeta, jumpzeta);
  double expzetastart = exp(zetastart);
  double alphastart = expzetastart/(1+expzetastart);
  
  // Declaration of some variable 
  double logalphazeta, logalpha2zeta;
  
  //Compute logalpha2
  List Astart(M), Aystart(M);
  
  arma::mat Gm = Gnorm[0];
  int Nm = N(0);
  arma::vec Aym = Ay[0];
  arma::vec ym = y[0];
  arma::vec Vmtheta  = Vtheta[0];
  
  arma::mat Astartm = arma::eye(Nm,Nm) - alphastart*Gm;
  arma::vec Aystartm = Astartm*ym;
  Astart[0] = Astartm; // save A*
  Aystart[0] = Aystartm; // save A*y
  double sumlogdetAstart = log(abs(det(Astartm))); // compute sumlogdetA
  
  logalpha2zeta = arma::dot(Aystartm,Aystartm) - arma::dot(Aym,Aym) +
    2*arma::dot(Vmtheta, Aym - Aystartm);
  
  for(int m(1); m<M; ++m) {
    arma::mat Gm = Gnorm[m];
    int Nm = N(m);
    arma::vec Aym = Ay[m];
    arma::vec ym = y[m];
    arma::vec Vmtheta  = Vtheta[m];
    
    Astartm = arma::eye(Nm,Nm) - alphastart*Gm;
    Aystartm = Astartm*ym;
    Astart[m] = Astartm;
    Aystart[m] = Aystartm;
    sumlogdetAstart += log(abs(det(Astartm)));
    
    logalpha2zeta += arma::dot(Aystartm,Aystartm) - arma::dot(Aym,Aym) +
      2*arma::dot(Vmtheta, Aym - Aystartm);
  }
  
  logalpha2zeta *= (-.5/sigma2);
  logalpha2zeta += sumlogdetAstart - sumlogdetA + 
    (0.5*invsigma2zeta)*(pow(zeta-zeta0,2) - pow(zetastart-zeta0,2));
  
  //Compute logalpha
  logalphazeta=min(Rcpp::NumericVector::create(0, logalpha2zeta));
  
  if(unif_rand()<exp(logalphazeta)){
    zeta = zetastart;
    alpha = alphastart;
    A = Astart;
    Ay = Aystart;
    sumlogdetA = sumlogdetAstart;
    zetaaccept +=1;     //Increase acceptance number to 1 
  }
}


// [[Rcpp::export]]
List peerMCMC (const List& y,
               const List& X,
               const List& Xone,
               List& Gnorm,
               const int& M,
               const IntegerVector& N,
               const int& kbeta,
               const int& kgamma,
               const List& prior,
               const arma::vec& theta0,
               const arma::mat& invsigmatheta,
               const double& zeta0,
               const double& invsigma2zeta,
               const double& a, 
               const double& b,
               const arma::vec parms0,
               int iteration,
               const double& target,
               const double& jumpmin,
               const double& jumpmax,
               const double& c,
               const bool& progress){
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  //N and Xshed, ListIndex
  List ListIndex(M);
  
  
  for(int m(0); m<M; ++m){
    // check prior 0 and 1;
    arma::mat priorm = prior[m];
    int Nm = N(m);
    arma::vec ListI(Nm);
    List ListIJ(Nm);
    int ci = 0;
    for(int i(0); i<Nm; ++i){
      int cj = 0;
      arma::vec ListJ(Nm);
      
      for(int j(0); j<Nm; ++j){
        if(priorm(i,j) !=0 && priorm(i,j) !=1){
          ListJ(cj) = j;
          cj += 1;
        }
      }
      if(cj > 0){
        ListI(ci) = i;
        ListJ = ListJ.head(cj);
        ListIJ[ci] = ListJ;
        ci += 1;
      }
    }
    
    ListI = ListI.head(ci);
    IntegerVector idx = seq_len(Nm);
    ListIJ = ListIJ[idx<=ci];
    ListIndex[m] = List::create(ListI,ListIJ);
  }
  
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
    sumlogdetA += log(abs(det(Am)));
  }
  
  
  
  //Save 
  arma::mat saveparms(kv+2,iteration);
  NumericVector parmscpp;
  // loop
  if (progress) {
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      
      // Update G
      updGnorm (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2, kbeta,kv,t);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      
      // update hyperparameters
      double jumpzetast = jumpzeta + (zetaaccept/t - target)/pow(t,c);
      if((jumpzetast > jumpmin) & (jumpzetast < jumpmax)){jumpzeta = jumpzetast;}
      
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      parmscpp             = wrap(parms);
      parmscpp.attr("dim") = R_NilValue;
      Rcpp::print(parmscpp);
      //Rcpp::Rcout<<"************************"<<std::endl;
      Rprintf("************************ \n");
    }
  } else {
    for(int t(0); t<iteration; ++t){
      // Update G
      updGnorm (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2, kbeta,kv,t);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      
      // update hyperparameters
      double jumpzetast = jumpzeta + (zetaaccept/t - target)/pow(t,c);
      if((jumpzetast > jumpmin) & (jumpzetast < jumpmax)){jumpzeta = jumpzetast;}
      
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
List peerMCMCblock (const List& y,
                    const List& X,
                    const List& Xone,
                    List& Gnorm,
                    const int& M,
                    const IntegerVector& N,
                    const int& kbeta,
                    const int& kgamma,
                    const List& prior,
                    const arma::vec& theta0,
                    const arma::mat& invsigmatheta,
                    const double& zeta0,
                    const double& invsigma2zeta,
                    const double& a, 
                    const double& b,
                    const arma::vec parms0,
                    int iteration,
                    const double& target,
                    const double& jumpmin,
                    const double& jumpmax,
                    const double& c, 
                    const int& nupmax,
                    const bool& progress){
  
  
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  //N and Xshed, ListIndex
  List ListIndex(M);
  
  
  for(int m(0); m<M; ++m){
    // check prior 0 and 1;
    arma::mat priorm = prior[m];
    int Nm = N(m);
    arma::vec ListI(Nm);
    List ListIJ(Nm);
    int ci = 0;
    for(int i(0); i<Nm; ++i){
      int cj = 0;
      arma::vec ListJ(Nm);
      
      for(int j(0); j<Nm; ++j){
        if(priorm(i,j) !=0 && priorm(i,j) !=1){
          ListJ(cj) = j;
          cj += 1;
        }
      }
      if(cj > 0){
        ListI(ci) = i;
        ListJ = ListJ.head(cj);
        ListIJ[ci] = ListJ;
        ci += 1;
      }
    }
    
    ListI = ListI.head(ci);
    IntegerVector idx = seq_len(Nm);
    ListIJ = ListIJ[idx<=ci];
    ListIndex[m] = List::create(ListI,ListIJ);
  }
  
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
    sumlogdetA += log(abs(det(Am)));
  }
  
  
  //Save 
  arma::mat saveparms(kv+2,iteration);
  NumericVector parmscpp;
  // loop
  if (progress) {
    for(int t(0); t<iteration; ++t){
      //std::cout<<"Iteration "<<t+1<<"/"<<iteration<<std::endl;
      Rprintf("Iteration %d/%d \n", t+1, iteration);
      // Update G
      updGnormblock (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2, kbeta,kv, nupmax);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      
      // update hyperparameters
      double jumpzetast = jumpzeta + (zetaaccept/t - target)/pow(t,c);
      if((jumpzetast > jumpmin) & (jumpzetast < jumpmax)){jumpzeta = jumpzetast;}
      
      arma::vec parms = arma::join_cols(theta,arma::ones(1)*alpha);
      parms = arma::join_cols(parms,arma::ones(1)*sigma2);
      
      saveparms.col(t) = parms;
      parmscpp             = wrap(parms);
      parmscpp.attr("dim") = R_NilValue;
      Rcpp::print(parmscpp);
      //Rcpp::Rcout<<"************************"<<std::endl;
      Rprintf("************************ \n");
    }
  } else {
    for(int t(0); t<iteration; ++t){
      // Update G
      updGnormblock (Gnorm, prior,ListIndex, N, M, y, A, Ay, Xb, Xgamma, alpha, sigma2, kbeta,kv, nupmax);
      for(int m(0); m<M; ++m){
        //Gnorm
        arma::mat Gm = Gnorm[m];
        
        // V
        arma::mat Xonem = Xone[m];
        arma::mat Xm = X[m];
        arma::mat Vm = arma::join_rows(Xonem,Gm*Xm);
        arma::vec Vmtheta = Vm*theta;
        V[m] = Vm;
      }
      
      // Update theta
      updtheta (theta, Vtheta, Xb, Xgamma, sigma2, kv, kbeta, Xone, X, Ay, V, invsigmathetatheta0,
                invsigmatheta, M);
      
      updsigma2 (sigma2, theta, a, b, theta0, invsigmatheta, Ay, Vtheta, sumN, M);
      
      updzeta (zeta, alpha, A, sumlogdetA, Ay, Gnorm, y, sigma2, Vtheta, jumpzeta,
               zetaaccept, zeta0, invsigma2zeta, N, M);
      
      // update hyperparameters
      double jumpzetast = jumpzeta + (zetaaccept/t - target)/pow(t,c);
      if((jumpzetast > jumpmin) & (jumpzetast < jumpmax)){jumpzeta = jumpzetast;}
      
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



class Peer: public Numer::MFuncGrad
{
private:
  List& listG;
  const Rcpp::IntegerVector& N;
  const int& M;
  const List& y;
  const List& X;
  const int& kz; //  number of parameters excepted sigma2
public:
  Peer(List& listG_, const Rcpp::IntegerVector& N_, const int& M_, const List& y_, const List& X_, const int& kz_) : 
  listG(listG_), N(N_), M(M_), y(y_), X(X_), kz(kz_){}
  
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
    
    List Gy = listG(0), GX = listG(1), EI = listG(2);
    
    arma::mat sZZ(kz,kz,arma::fill::zeros);
    arma::vec sZAY(kz,arma::fill::zeros);
    arma::vec sZGY(kz,arma::fill::zeros);
    
    double sumN = sum(N);
    //compute beta as [beta,gamma]
    
    for(int m(0); m<M; ++m){
      arma::mat Xm = X[m],  GXm = GX[m];
      arma::vec Gym = Gy[m], ym=y[m];
      
      arma::mat Zm = arma::join_rows(Xm,GXm);
      
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
      arma::mat Xm = X[m],  GXm = GX[m];
      arma::vec Gym = Gy[m], ym=y[m];
      
      arma::mat Zm = arma::join_rows(Xm,GXm);
      
      sA += (arma::trans(alpha*Gym - ym)*Gym);
      sB += Gym.t()*Zm;
      sC += Zm.t()*Zm;
      sD += ym.t()*Zm;
      
      arma::vec em = ym - alpha*Gym - Zm*beta;
      
      ee += dot(em,em);
    }
    
    sigma2 = ee/sumN;
    dsigma2 = arma::conv_to<double>::from(sA + sB*(beta + alpha*dbeta) + (beta.t()*sC - sD)*dbeta);
    
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
                    const List& X,
                    const int& kv) {
  // Negative log likelihood
  Peer f(listG, N, M, y, X, kv);
  
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
  List Gy(M), GX(M), EIOGGX(M);
  
  for(int m(0); m<M; ++m){
    arma::mat Gm = Gnorm(m);
    arma::mat Xm = Xone(m);
    arma::mat Xmshed = X(m);
    arma::vec ym = y(m);
    
    // GY
    Gy(m) = Gm*ym;
    
    // GX
    GX(m) = Gm*Xmshed;
    
    // Eigen
    arma::cx_vec EIGENm;
    arma::eig_gen(EIGENm, Gm);
    EIOGGX(m) = EIGENm;
  }
  
  listG(0) = Gy;
  listG(1) = GX;
  listG(2) = EIOGGX;
  
  List MLestim=estimatepeerML(listG, N, M, y, Xone, kz);
  
  double alphaML = MLestim(0);
  arma::vec thetaML = MLestim(1);
  double seML = MLestim(2);
  
  theta.head(kz) = thetaML;
  theta(kz)=alphaML;
  theta(kz+1)=seML;
  
  return theta;
}


