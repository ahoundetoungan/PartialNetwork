// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppEigen.h>
#include "vMF.h"

using namespace Rcpp;
//using namespace arma;
using namespace std;
//using namespace Eigen;
//using namespace Numer;


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
    const arma::Col<int>& index_col,
    const arma::rowvec& newval,
    const int& nupdate
){
  
  for (int j(0); j < nupdate; ++ j) {
    Gmnorm(i, index_col(j)) = newval(j);
    prior_blockl *= (priorm(i, index_col(j))*newval(j) + (1-priorm(i, index_col(j)))*(1-newval(j)));
  }
}

// compute the entry of the matrix to update given the prior density  and the part of the network which is observed
// [[Rcpp::export]]
List fListIndex(List& prior,
                List& G0obs,
                const int& M,
                const IntegerVector& N) {
  List ListIndex(M);
  for(int m(0); m<M; ++m){
    arma::mat G0obsm = G0obs[m];
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
        if(priorm(i,j) !=0 && priorm(i,j) !=1 && G0obsm(i, j) <= 0){
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
  return ListIndex;
}

// This function normalizes network
// [[Rcpp::export]]
List fGnormalise(List& u, const double& M) {
  List out(M);
  
  for(int m(0); m < M; ++m) {
    arma::mat um = u[m];
    um           = arma::normalise(um, 1, 1);
    out[m]       = um;
  }
  return out;
}


// Create a list a square matrixes from a given vector and sizes.
// The size of the m-th matrice in the list is N[m] with zeros on the diagonal
// The elements in the generated matrix are placed column-wise (ie. the first column is filled up before filling the second column)
// [[Rcpp::export]]
List frVtoM(const Eigen::VectorXd& u,
            const Rcpp::IntegerVector& N,
            const double& M) {
  List out(M);
  
  int r                                = 0;
  int n;
  
  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    
    n                                  = Nm - 1;
    
    Eigen::MatrixXd outm(Eigen::MatrixXd::Zero(Nm, Nm));
    outm.block(1, 0, n, 1)             = u.segment(r, n);
    
    r                                 += n;
    for(int i(1); i < n; ++i) {
      outm.block(0, i, i, 1)          = u.segment(r, i);
      outm.block(i + 1, i, n - i, 1)  = u.segment(r + i, n - i);
      r                              += n;
    }
    
    outm.block(0, n, n, 1)            = u.segment(r, n);
    r                                += n;
    
    out[m]                            = outm;
  }
  return out;
}

// does the same thing but the entry matrix is armadillo
List frVtoMarma(const arma::vec& u,
                const Rcpp::IntegerVector& N,
                const double& M) {
  List out(M);
  
  int r                              = 0;
  int n;
  
  for(int m(0); m < M; ++m) {
    int Nm                           = N(m);
    
    n                                = Nm - 1;
    
    arma::mat outm(Nm, Nm, arma::fill::zeros);
    outm.submat(1, 0, n, 0)          = u.subvec(r, n + r - 1);
    
    r                               += n;
    for(int i(1); i < n; ++i) {
      outm.submat(0, i, i - 1, i)    = u.subvec(r, r + i - 1);
      outm.submat(i + 1, i, n, i)    = u.subvec(r + i, n + r - 1);
      r                             += n;
    }
    
    outm.submat(0, n, n - 1, n)      = u.subvec(r, r + n - 1);
    r                               += n;
    
    out[m]                           = outm;
  }
  return out;
}

// Same function but the returned matrix are normalized
// [[Rcpp::export]]
List frVtoMnorm(const arma::vec& u,
                const IntegerVector& N,
                const double& M) {
  List out(M);
  
  int r2                               = -1;
  int r1;
  
  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    
    r2                                += Nm - 1;
    r1                                 = r2 - Nm + 2;
    
    arma::mat outm(Nm, Nm, arma::fill::zeros);
    outm.submat(1, 0, Nm - 1, 0)       = u.subvec(r1, r2);
    
    for(int i(1); i < (Nm - 1); ++i) {
      r2                              += Nm - 1;
      r1                               = r2 - Nm + 2;
      outm.submat(0, i, i - 1, i)      = u.subvec(r1, r1 + i - 1);
      outm.submat(i + 1, i, Nm - 1, i) = u.subvec(r1 + i, r2);
    }
    
    r2                                += Nm - 1;
    r1                                 = r2 - Nm + 2;
    outm.submat(0, Nm - 1, Nm - 2, Nm - 1) = u.subvec(r1, r2);
    
    outm                                   = arma::normalise(outm, 1, 1);
    
    out[m]                                 = outm;
  }
  return out;
}


// Create a vector from a given list a square matrixes
// The size of the length of the vector is the sum(N), where N is the vector of matrice sizes
// The elements in the generated vector are taken from column-wise (ie. the first column is filled up before filling the second column)
// and from the first matrix of the list to the last matrix of the list.
// [[Rcpp::export]]
Eigen::VectorXd frMtoV(List& u,
                       const Rcpp::IntegerVector& N,
                       const double& M) {
  int sN                               = sum(N*N - N);
  Eigen::VectorXd out(sN);
  
  int r                                = 0;
  int n;
  
  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    Eigen::MatrixXd um                 = u[m];
    //um                                 = um.array().ceil();
    
    n                                  = Nm - 1;
    
    out.segment(r, n)                  = um.block(1, 0, n, 1);    
    r                                 += n;
    for(int i(1); i < n; ++i) {
      out.segment(r, i)                = um.block(0, i, i, 1);
      out.segment(r + i, n - i)        = um.block(i + 1, i, n - i, 1);
      r                               += n;
    }
    
    out.segment(r, n)                  = um.block(0, n, n, 1);
    r                                 += n;
  }
  return out;
}

// same function but the matrixes are ceiled first
// [[Rcpp::export]]
Eigen::VectorXd frMceiltoV(List& u,
                           const Rcpp::IntegerVector& N,
                           const double& M) {
  int sN                               = sum(N*N - N);
  Eigen::VectorXd out(sN);
  
  int r                                = 0;
  int n;
  
  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    Eigen::MatrixXd um                 = u[m];
    um                                 = um.array().ceil();
    
    n                                  = Nm - 1;
    
    out.segment(r, n)                  = um.block(1, 0, n, 1);    
    r                                 += n;
    for(int i(1); i < n; ++i) {
      out.segment(r, i)                = um.block(0, i, i, 1);
      out.segment(r + i, n - i)        = um.block(i + 1, i, n - i, 1);
      r                               += n;
    }
    
    out.segment(r, n)                  = um.block(0, n, n, 1);
    r                                 += n;
  }
  return out;
}

// This function set the jumping scales in min and max
// vector list version
void fsetjump_vl(arma::vec& jump, const double& jumpmin,
                const double& jumpmax, const int& M,
                List& jumpl, const List& Vparms) {
  arma::uvec acc_updmin = arma::find(jump < jumpmin);
  jump.rows(acc_updmin) = arma::ones(acc_updmin.n_elem)*jumpmin;
  arma::uvec acc_updmax = arma::find(jump > jumpmax);
  jump.rows(acc_updmax) = arma::ones(acc_updmax.n_elem)*jumpmax;
  for(int m(0); m < M; ++ m) {
    const arma::mat Vrhom    = Vparms[m];
    jumpl[m]                 = jump(m)*jump(m)*Vrhom;
  }
}
// double matrix version
void fsetjump_dm(double& jump, const double& jumpmin,
                const double& jumpmax, arma::mat& jumpl, const arma::mat& Vparms) {
  if(jump < jumpmin) jump = jumpmin;
  if(jump > jumpmax) jump = jumpmax;
  jumpl                   = jump*jump*Vparms;
}
// vec version
void fsetjump_v(arma::vec& jump, const double& jumpmin,
                const double& jumpmax) {
  arma::uvec acc_updmin = arma::find(jump < jumpmin);
  jump.rows(acc_updmin) = arma::ones(acc_updmin.n_elem)*jumpmin;
  arma::uvec acc_updmax = arma::find(jump > jumpmax);
  jump.rows(acc_updmax) = arma::ones(acc_updmax.n_elem)*jumpmax;
}
// row version
void fsetjump_r(arma::rowvec& jump, const double& jumpmin,
                const double& jumpmax) {
  arma::uvec acc_updmin = arma::find(jump < jumpmin);
  jump.cols(acc_updmin) = arma::ones<arma::rowvec>(acc_updmin.n_elem)*jumpmin;
  arma::uvec acc_updmax = arma::find(jump > jumpmax);
  jump.cols(acc_updmax) = arma::ones<arma::rowvec>(acc_updmax.n_elem)*jumpmax;
}
// row double
void fsetjump_d(double& jump, const double& jumpmin,
                const double& jumpmax) {
  if(jump < jumpmin) jump = jumpmin;
  if(jump > jumpmax) jump = jumpmax;
}

////// Estimate graph Breza et al. 

// Compute distance matrix between people
// Distance between individual described by catagorial variables. We use
// J. Pagès [2004] "Analyse Factorielle des Données mixtes", Revue de 
// Statistique Appliquée, tome 52, N°4, P. 93-111

// after we can build the index and the weights matrix to compute nu and  
// z for people without ARD. The matrix dimension is N2*m
// The entry [m,i] is m-th neighbor for i in index matrix (col(i) are the m neighbors of i)
// and its weights in weight matrix defined by 1/(1+d[j,i]), where row(i) are the weight of i

void cneighbor(const double& N1, const double& N2, const double& N,
               const arma::mat& Xard, const arma::mat& Xnonard, 
               const int& m, arma::umat &neighbor, arma::mat& weight){
  
  //Normalize disjuntif tables following Pagès [2004]
  arma::mat traitnorm = arma::normalise(arma::join_cols(Xard,Xnonard),2,0)*sqrt(N);
  
  //variable for saving distance of ARD j to ARD i
  arma::vec disti(N1);
  arma::mat dist(m, N2);
  //other variables
  arma::uvec temp, index;         
  
  //Compute distance
  for(int i(0);i<N2;++i){
    for(int j(0);j<N1;++j){
      // j is ARD and i nonARD
      disti(j) = arma::norm(traitnorm.row(j)-traitnorm.row(N1+i));
    }
    temp            = arma::sort_index(disti);        //Returns index for distance from smallest 
    index           = temp.head(m);                   //the m-th nearest index 
    neighbor.col(i) = index;
    dist.col(i)     = disti.elem(index);
  }
  
  dist              = 1.0/(1.0 + dist.t());
  weight            = arma::normalise(dist, 1, 1);
}

// This function add non ARD objects to ARD objects
// ARD objects are zm, num and dm
void frhononARD(arma::mat& zm, arma::vec& num, arma::vec& dm, const double& logCpzeta,  const int& N1m, const int& N2m,
                const int& Nm,  const int& P, const arma::umat& neighbor, const arma::mat& weight, 
                const arma::uvec& iARD, const arma::uvec& inonARD) {
  arma::mat zall(Nm, P), znonARD(N2m, P);
  zall.rows(iARD)     = zm;
  arma::vec nuall(Nm), dall(Nm), nunonARD(N2m);
  nuall.elem(iARD)    = num;
  dall.elem(iARD)     = dm;
  
  arma::uvec neighborj;
  arma::rowvec weightj;
  //compute nu and z for non ARD
  for(int j(0); j < N2m; ++j){
    neighborj         = neighbor.col(j);
    weightj           = weight.row(j);
    nunonARD.row(j)   = weightj*(num.elem(neighborj));
    znonARD.row(j)    = weightj*(zm.rows(neighborj));
  }
  znonARD             = arma::normalise(znonARD, 2, 1);
  zall.rows(inonARD)  = znonARD;
  nuall.rows(inonARD) = nunonARD;
  dall.rows(inonARD)  = (Nm*1.0/N1m)*arma::exp(nunonARD)*arma::accu(exp(num))/exp(logCpzeta);
  
  // compute Probabilities
  zm                  = zall;
  num                 = nuall;
  dm                  = dall;
}

// This function returns zeta, z, nu, d, ie  (initial values for sarARD)
// [[Rcpp::export]]
List flspacerho1(const double& T, const double& P, const arma::cube& z, const arma::mat& d,
                      const arma::vec& zeta, const unsigned int& N1,  const unsigned int& Metrostart){
  // Number of people with ARD
  double N = N1, zetat, logCpzetat, ngraph = T-Metrostart+1;      
  arma::rowvec nuARDt(N), ds(N, arma::fill::zeros), dt;
  arma::mat zt;
  
  //export
  arma::mat znut(N, P + 1);
  arma::mat znu(ngraph, N*(P + 1) + 1);
  znu.col(0)        = zeta.tail(ngraph);
  
  //loop 
  for(int t(0);t<ngraph;++t){  
    zetat           = zeta(t+Metrostart);      // extract zeta for iteration t+Metrostart
    zt              = z.slice(t+Metrostart);   // extract z for iteration t+Metrostart
    dt              = d.row(t+Metrostart);     // extract d for iteration t+Metrostart
    ds             += dt;
    logCpzetat      = logCpvMFcpp(P,zetat);     
    //compute nu for ARD
    nuARDt          = log(dt) + 0.5*logCpzetat - 0.5*log(arma::accu(dt)) ;
    
    //save
    znut                           = arma::join_rows(zt, nuARDt.t());
    znut.reshape(1, N*(P + 1));
    znu.submat(t, 1, t, N*(P + 1)) = znut;
  }
  ds   /= ngraph;
  return List::create(Named("rho")    = znu,
                      Named("degree") = ds);
}
// [[Rcpp::export]]
List flspacerho2(const double& T, const double& P, const arma::cube& z, const arma::mat& d,
                      const arma::vec& zeta, const arma::mat& Xard, const arma::mat& Xnonard, 
                      const unsigned int& N1, const unsigned int& N2, const unsigned int& M, 
                      const unsigned int& Metrostart){
  // Number of people with ARD
  double N = N1 + N2, zetat, logCpzetat,  ngraph = T - Metrostart + 1;
  arma::rowvec nut, dt, ds(N1, arma::fill::zeros);
  arma::mat zt;
  
  //export
  arma::mat znu(ngraph, N1*(P + 1) + 1);
  znu.col(0)        = zeta.tail(ngraph);
  
  //compute neighbor and weight
  arma::umat neighbor(M, N2);            
  arma::mat weight; // will the transposed in the function cneighbor and output will be (N2, M) 
  cneighbor(N1, N2, N, Xard, Xnonard, M, neighbor, weight);
  
  //Necessary variables
  arma::uvec neighborj;
  arma::rowvec weightj;
  
  //loop 
  for(int t(0);t<ngraph;++t){  // thus convergence for t from Metrostart to T+1
    zetat           = zeta(t+Metrostart);     
    zt              = z.slice(t+Metrostart);       
    dt              = d.row(t+Metrostart); 
    ds             += dt;
    logCpzetat      = logCpvMFcpp(P,zetat);     
    //compute nu for ARD
    nut             = log(dt) + 0.5*logCpzetat + 0.5*log(N1*1.0/N) - 0.5*log(arma::accu(dt)) ;
    
    //save
    arma::mat znut  = arma::join_rows(zt, nut.t());
    znut.reshape(1, N1*(P + 1));
    znu.submat(t, 1, t, N1*(P + 1)) = znut;
  }
  ds   /= ngraph;
  return List::create(Named("rho")      = znu,
                      Named("degree")   = ds,
                      Named("neighbor") = neighbor,
                      Named("weights")   = weight);
}

// Given ZNU this function compute the probability of Gij = 1 using the latent space model
// [[Rcpp::export]]
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
                   const arma::uvec& inonARDm) {
  
  if(N2m > 0) {
    frhononARD(zm, num, dm, logCpzetam, N1m, N2m, Nm,  Pm, neighborm, weightm, iARDm, inonARDm);
  }
  zm              = arma::normalise(zm, 2, 1);
  arma::mat numat = arma::repmat(num, 1, Nm);
  arma::mat prob  = arma::exp(zetam*zm*zm.t() + numat + numat.t());
  
  prob           *= ((sum(dm))/arma::accu(prob));
  arma::uvec tmp  = arma::find(prob > 1);
  prob.elem(tmp)  = arma::ones(tmp.n_elem);
  prob.diag()     = arma::zeros(Nm);
  return prob;
}

// This function creates a list of armadillo matrix
// M is the length of the list
// Krho is a vector of the number of columns in each matrix
List createlistmat(const int& M, const arma::vec& Krho) {
  List out(M);
  for (int m(0); m < M; ++m) {
    arma::mat outm(0, Krho(m));
    out[m] = outm;
  }
  return out;
}

// This function add rows to matrix in a list
// M is the length of the list
// listmat is the current list mat
// rho is a list of the rows to add to each matrix
void addtolistmat(const int& M, List& listmat, List& Rho) {
  for (int m(0); m < M; ++m) {
    arma::vec Rhom     = Rho[m];
    arma::mat listmatm = listmat[m];
    listmatm           = arma::join_cols(listmatm, Rhom.t());
    listmat[m]         = listmatm;
  }
}
