// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace Eigen;
using namespace Numer;

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

//[[Rcpp::export]]
Eigen::MatrixXd invmodij (const Eigen::MatrixXd& invM,
                          const unsigned int& i,
                          const unsigned int& j,
                          const double& eps_) {
  return invM - eps_*invM.col(i)*invM.row(j)/(1. + eps_*invM(j,i));
}

//[[Rcpp::export]]
Eigen::MatrixXd invmodijk (const Eigen::MatrixXd& invM,
                           const unsigned int& i,
                           const unsigned int& j,
                           const unsigned int& k,
                           const double& eps_1,
                           const double& eps_2) {
  Eigen::MatrixXd tmp = invM - eps_1*invM.col(i)*invM.row(j)/(1 + eps_1*invM(j,i));
  return tmp - eps_2*tmp.col(i)*tmp.row(k)/(1 + eps_2*tmp(k,i));
}

//[[Rcpp::export]]
double detmodij (const double& detM,
                 const Eigen::MatrixXd& invM,
                 const unsigned int& i,
                 const unsigned int& j,
                 const double& eps_) {
  return (1 + eps_*invM(j,i))*detM;
}

//[[Rcpp::export]]
double detmodijk (const double& detM,
                  const Eigen::MatrixXd& invM,
                  const unsigned int& i,
                  const unsigned int& j,
                  const unsigned int& k,
                  const double& eps_1,
                  const double& eps_2) {
  return (1 + eps_2*invM(k,i) - eps_1*eps_2*invM(k,i)*invM(j,i))*(1 + eps_1*invM(j,i))*detM;
}

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
              double& zeta0,
              double& invsigma2zeta,
              const Rcpp::IntegerVector N,
              const double M) {
  
  NumericVector zetatemp = rnorm(1, zeta, jumpzeta);
  
  double zetastart = zetatemp(0);
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
NumericMatrix peerMCMC (const List& y,
                        const List& X,
                        const arma::vec parms0,
                        const List& hyperparms,
                        const int& iteration = 1,
                        const bool& intercept = true,
                        const double& target = 0.44,
                        const double& jumpmin = 1e-12,
                        const double& jumpmax = 100,
                        const double& c = 0.6){
  
  clock_t begintime=clock();  
  List Xone = clone(X); // Clone X
  
  // hyper parameters
  const List prior = hyperparms(0);
  arma::vec theta0 = hyperparms(1);
  const arma::mat invsigmatheta = hyperparms(2);
  double zeta0 = hyperparms(3);
  double invsigma2zeta = hyperparms(4);
  double a = hyperparms(5);
  double b = hyperparms(6);
  const arma::vec& invsigmathetatheta0 = invsigmatheta*theta0;
  
  //M
  int M = y.length();
  
  //N and Xshed, ListIndex
  Rcpp::IntegerVector N(M);
  List ListIndex(M);
  
  
  for(int m(0); m<M; ++m){
    arma::mat Xm = Xone(m);
    N(m) = Xm.n_rows;
    
    int Nm = N(m);
    if(intercept){
      Xm = arma::join_rows(arma::ones(Nm),Xm);
      Xone(m) = Xm;
    }
    
    // check prior 0 and 1;
    arma::mat priorm = prior[m];
    
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
  
  //kz
  int kbeta, kgamma;
  arma::mat X0 = Xone(0);
  kbeta = X0.n_cols;
  
  kgamma = kbeta;
  
  if(intercept){
    kgamma -= 1;
  }
  
  int kv = kbeta + kgamma; // number of exogenous variables
  
  //initialize parameters
  arma::vec theta = parms0.head(kv);
  double sigma2 = parms0(kv+1);
  double alpha = parms0(kv);
  double zeta = log(alpha/(1-alpha));
  
  // Other parameters
  double jumpzeta = 1;
  double zetaaccept = 0;
  List Gnorm(M), V(M), Vtheta(M), Xb(M), Xgamma(M), A(M), Ay(M);
  double sumlogdetA = 0.0;
  for(int m(0); m<M; ++m){
    //Gnorm
    int Nm = N(m);
    arma::mat priorm = prior[m];
    mat matunif(Nm,Nm,fill::randu);
    arma::mat Gm = arma::normalise(arma::conv_to<arma::mat>::from(matunif < priorm),1,1);
    Gnorm[m] = Gm;
    
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

  // loop
  for(int t(0); t<iteration; ++t){
    cout<<"Iteration "<<t+1<<"/"<<iteration<<endl;
    
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
    cout<<parms.t()<<endl;
  }
  
  Rcpp::Environment base("package:base"); 
  
  // Make function callable from C++
  Rcpp::Function cCPP = base["c"];  
  Rcpp::Function paste0CPP = base["paste0"];  
  
  StringVector Names = cCPP(paste0CPP("X",seq_len(kgamma)),paste0CPP("GX",seq_len(kgamma)),"alpha","sd.");
  if(intercept){
    Names = cCPP("Intercept",Names);
  }
  
  NumericMatrix output = wrap(saveparms.t());
  colnames(output) = Names;
  
  //Print the processing time
  double timer=((double)(clock()-begintime))/CLOCKS_PER_SEC;//Time in seconds
  long nhours=floor(timer/3600);                            //Hours
  long nminutes=((long)floor(timer/60))%60;                 //Minutes
  double nseconds=timer-3600*nhours-60*nminutes;            //Seconds
  cout<<"Elapsed time   : "<<nhours<<" HH "<<nminutes<<" mm "<<round(nseconds)<<" ss "<<endl;
  cout<<""<<endl;
  cout<<"zeta average acceptance rate: "<<zetaaccept/iteration<<endl;
  
  return output;
}


class Peer: public MFuncGrad
          {
private:
  List& listG;
  const Rcpp::IntegerVector& N;
  const int& M;
  const List& y;
  const List& X;
  const int& kz; //  number of parameters expected sigma2
public:
  Peer(List& listG_, const Rcpp::IntegerVector& N_, const int& M_, const List& y_, const List& X_, const int& kz_) : 
  listG(listG_), N(N_), M(M_), y(y_), X(X_), kz(kz_){}
  
  double alpha;
  arma::vec beta;
  double sigma2;
  arma::vec dbeta;
  double dsigma2;
  double Grad;
  
  
  double f_grad(Constvec& parms, Refvec grad)
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
  
  //Temporary
  cout<<"Objective function: "<< -fopt<<endl;
  //End Temporary
  
  return Rcpp::List::create(
    Rcpp::Named("alpha") = f.alpha,
    Rcpp::Named("beta") = f.beta,
    Rcpp::Named("sigma2") = f.sigma2,
    Rcpp::Named("status") = status);
}


// Starting point
// [[Rcpp::export]] 
NumericVector sartpoint (const List& prior,
                         const List& y,
                         const List& X,
                         const bool& intercept = true) {
  
  List Xone = clone(X);
  
  //M
  int M = y.length();
  
  //N and Xshed
  Rcpp::IntegerVector N(M);
  
  for(int m(0); m<M; ++m){
    arma::mat Xm = Xone(m);
    N(m) = Xm.n_rows;
    
    if(intercept){
      Xm = arma::join_rows(arma::ones(N(m)),Xm);
      Xone(m) = Xm;
    }
  }
  
  //kz
  int kbeta, kgamma;
  arma::mat X0 = Xone(0);
  kbeta = X0.n_cols;
  
  kgamma = kbeta;
  
  if(intercept){
    kgamma -= 1;
  }
  
  int kz = kbeta + kgamma; // number of exogenous variables
  
  // Initialize G
  List Gnorm(M);
  for(int m(0); m<M; ++m){
    arma::mat priorm = prior(m);
    mat matunif(N(m),N(m),fill::randu);
    arma::mat G = arma::conv_to<arma::mat>::from(matunif < priorm);
    Gnorm(m) = arma::normalise(G,1,1);
  }
  
  
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

  Rcpp::Environment base("package:base"); 
  
  // Make function callable from C++
  Rcpp::Function cCPP = base["c"];  
  Rcpp::Function paste0CPP = base["paste0"];  
  
  StringVector Names = cCPP(paste0CPP("X",seq_len(kgamma)),paste0CPP("GX",seq_len(kgamma)),"alpha","sd.");
  if(intercept){
    Names = cCPP("Intercept",Names);
  }
  
  NumericVector output = wrap(theta);
  
  rownames(output) = Names;
  
  return output;
}


