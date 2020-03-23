// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppEigen.h>

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
    const arma::vec& index_col,
    const arma::rowvec& newval,
    const int& nupdate
){
  
  for (int j(0); j < nupdate; ++ j) {
    Gmnorm(i, index_col(j)) = newval(j);
    prior_blockl *= (priorm(i, index_col(j))*newval(j) + (1-priorm(i, index_col(j)))*(1-newval(j)));
  }
}