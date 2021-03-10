// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <progress.hpp>
#include <progress_bar.hpp>
#include "tools.h"
#include "vMF.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

/*
 *  Estimate graph parameters Mccormick and Zheng (2015)
 *  
 *  To sample from von Mises Fisher Distribution
 *  We follow Wood, A. T. (1994). Simulation of the von Mises Fisher distribution.
 *  Communications in statistics-simulation and computation, 23(1), 157-164.
 *  
 *  In a first step we built a function rw to compute W; that is, 
 *  perform steps 0, 1 et 2.
 *  
 *  The function needs the sample size "size", lambda as intensity parameter,
 *  d as dimension p-1 et a Vector W (void) of dimension "size" to save the
 *  values corresponding to each sampling.
 *  
 *  In the Metropolis steps, one zi is draw for each i and for each iteration
 *  Thus we build a fast version of previous functions for size = 1.
 *  That allows to avoid the loop over n and is faster.
 *  This function is not exported because it is used only in the C++ code.
 */

//****** Some useful functions
// Compute matrix Cp with Cpik=Cp(norm(zeta*zi+etak*vk))
arma::mat computelogCp(const double& n, const double& K, const double& p,
                       const arma::mat& z, const arma::mat& v, const arma::rowvec& eta,
                       const double& zeta){
  
  arma::rowvec etaeta = eta%eta, temp;
  arma::mat zvt       = z*v.t();
  double zetazeta     = zeta*zeta;
  
  arma::mat Output(n,K);
  
  for(int i(0); i < n; ++ i){
    temp = sqrt(zetazeta + etaeta + 2*zeta*eta%zvt.row(i));
    for(int k(0); k < K; ++ k){
      Output(i,k) = logCpvMFcpp(p, temp(k));
    }
  }
  
  return Output;
}


arma::mat computelogL(const int& n, const int& K, const arma::vec& logd, const arma::rowvec& logb,
                      const arma::rowvec& logCpeta,  const double& logCpzeta, const arma::mat& logCp){
  arma::mat llambda(n, K);
  arma::rowvec tmp = logb + logCpeta + logCpzeta;
  
  for (int k(0); k < K; ++ k) {
    llambda.col(k) = tmp(k) + logd;
  }
  llambda -= logCp;
  return llambda;
}



// Compute the matrix L with Lik=lambdaik
arma::mat loglikelihood(const int& n, const int& K, const arma::vec& logd, const arma::rowvec& logb,
                        const arma::rowvec& logCpeta,  const double& logCpzeta, const arma::mat& logCp,
                        const arma::mat& Y){
  arma::mat llambda = computelogL(n, K,  logd, logb, logCpeta, logCpzeta, logCp);
  arma::mat lambda  = exp(llambda);
  
  return -lambda + Y%llambda;
}


// Six functions to update the six parameters : z, v, d, b, eta and zeta
// The argyments are
//n: sample size
//K: the number of groups
//p: the hypersphere dimension
//Y: the matrix of ARD
//z, v, d, b, eta, zeta: current values of the parameters
//jump[z, d, b, eta, zeta]: jump scaling of th eparameters
//[z, d, b, eta, zeta]accept: to count if the candidate value has been accept
//Cpeta: a vector of Cp(eta(k))
//Cpzeta: a scalar Cp(zeta)
//phieta: a vector of phi(eta(k)/jumpeta)
//phizeta: ascalar phi(zeta/jumpzeta)
//Cp: the matrix Cp
//L: the matrix L

//This step follows the section 3.3 in the technical detail
//The functions do not return a value. They take just the current values
//And change their values if the proposals are accepted

//Update z

void zupdate(const double& n, const double& K, const double& p,const arma::mat& Y,
             arma::mat& z,  const arma::mat& v, const arma::vec& logd, const arma::rowvec& logb, 
             const arma::rowvec& eta, const double& zeta, const arma::vec& jumpz, arma::vec& zaccept,
             const arma::rowvec& logCpeta, double& logCpzeta, arma::mat& logCp, arma::mat& loglik){
  //Draw from the proposal
  arma::mat zstart(n,p);
  for(int i(0);i<n;++i){
    zstart.row(i)=rvMFone(p,jumpz(i)*(z.row(i)).t()) ; 
  }
  //Step 1: Compute Cpistart and Listart
  arma::mat logCpstart  = computelogCp(n, K, p, zstart, v, eta, zeta);
  arma::mat loglikstart = loglikelihood(n, K, logd, logb, logCpeta, logCpzeta, logCpstart, Y);
  //Compute logalpha
  arma::vec logalpha    = min(arma::join_rows(arma::zeros(n), sum(loglikstart - loglik, 1)), 1);
  arma::vec r_unif(n, arma::fill::randu);
  arma::uvec selupd     = arma::find(r_unif < exp(logalpha));
  // update
  z.rows(selupd)        = zstart.rows(selupd);  
  logCp.rows(selupd)    = logCpstart.rows(selupd);    
  loglik.rows(selupd)   = loglikstart.rows(selupd);
  zaccept.elem(selupd) += arma::ones(selupd.n_elem); //Increase acceptance number to 1 
}


//Update v

void vupdate(const double& K, const double& p, List& indexg, const arma::mat& z, arma::mat& v,
             const arma::rowvec& eta,  const arma::uvec& fixv, const arma::mat& vf){
  
  //To save the sum of zi with i in Gk
  arma::rowvec sumzk;
  
  //Loop over k
  for(int k(0);k<K;++k){
    //Compute the sum of zi with i in Gk
    sumzk= sum(z.rows(as<arma::uvec>(indexg(k))),0);
    //Draw vk
    v.row(k)=rvMFone(p,eta(k)*sumzk.t());
  }
  
  
  //Rotate back for fixed column
  //Extract values to rotate. We call it vs (see section 3.3)
  arma::mat vs=v.rows(fixv-1);
  //Compute A
  arma::mat A=vs.t()*vf*vf.t()*vs;
  //compute A^[1/2]. We call it A2
  arma::vec Aeval;
  arma::mat Aevec;
  eig_sym(Aeval,Aevec,A);
  arma::mat A2=Aevec*diagmat(sqrt(Aeval))*(Aevec).t();
  //Ratation back
  v.rows(fixv-1)=vf;//vs*inv(A2)*vs.t()*vf;
}


//Update d

void dupdate(const int& n, const int& K, const arma::mat& Y, const arma::mat& z, const arma::mat& v, arma::vec& logd, 
             const arma::rowvec& logb, const arma::rowvec& eta, const double& zeta, const double& mud,
             const double& sigmad, const arma::vec& jumpd, arma::vec& daccept, const arma::rowvec& logCpeta, 
             const double& logCpzeta, const arma::mat& logCp, const arma::vec& SumYrow){
  
  // Declaration of some variable 
  arma::mat tmp = logCpzeta - logCp;
  arma::rowvec tmprv = logb + logCpeta;
  for (int k(0); k < K; ++ k) {
    tmp.col(k) += tmprv(k);
  }
  // Draw from the proposal
  arma::vec logdstart  = logd + jumpd%arma::randn<arma::vec>(n); 
  //Compute logalpha
  arma::vec logalpha   = 0.5*(pow((logd - mud)/sigmad, 2) - pow((logdstart - mud)/sigmad, 2))  +
    (exp(logd)-exp(logdstart))%sum(exp(tmp),1) + 
    SumYrow%(logdstart - logd);
  arma::vec r_unif(n, arma::fill::randu);
  arma::uvec selupd    = arma::find(r_unif < exp(logalpha));
  // update
  logd.elem(selupd)     = logdstart.elem(selupd);  
  daccept.elem(selupd) += arma::ones(selupd.n_elem); //Increase acceptance number to 1 
}



//Update b
void bupdate(const int& n, const double& K,const arma::mat& Y, const arma::mat& z, const arma::mat& v,
             const arma::vec& logd,  arma::rowvec& logb, const arma::rowvec& eta, const double& zeta,
             const double& mub, const double& sigmab,  const arma::rowvec& jumpb, 
             arma::rowvec& baccept, const arma::rowvec& logCpeta, const double& logCpzeta,
             const arma::mat& logCp, const arma::rowvec& SumYcol){
  
  // Declaration of some variable 
  arma::mat tmp  = logCpzeta - logCp;
  for (int i(0); i < n; ++ i) {
    tmp.row(i)  += logd(i);
  }
  
  // Draw from the proposal
  arma::rowvec logbstart = logb + jumpb%arma::randn<arma::rowvec>(K); 
  //Compute logalpha
  arma::rowvec logalpha  = (0.5/(sigmab*sigmab))*(logb - logbstart)%(logb - logbstart - 2*mub) +
    (exp(logb) - exp(logbstart))%exp(logCpeta)%sum(exp(tmp),0) + 
    SumYcol%(logbstart - logb);
  arma::vec r_unif(K, arma::fill::randu);
  arma::uvec selupd    = arma::find(r_unif < exp(logalpha.t()));
  // update
  logb.elem(selupd)     = logbstart.elem(selupd);  
  baccept.elem(selupd) += arma::ones<arma::rowvec>(selupd.n_elem); //Increase acceptance number to 1 
}

//Update eta

void etaupdate(const double& n, const double& K, const double& p, const arma::mat& Y,
               arma::mat& z, arma::mat& v, arma::vec& logd, arma::rowvec& logb, arma::rowvec& eta, double& zeta,
               const double& alphaeta, const double& betaeta, 
               const arma::rowvec& jumpeta, arma::rowvec& etaaccept, arma::rowvec& logCpeta,
               double& logCpzeta, arma::mat& logCp, arma::mat& loglik){
  
  //Variable of the candidate value
  arma::rowvec etastart(K); 
  //Declaration of some variable 
  arma::rowvec logCpetastart(K);
  
  //Loop over k
  for(int k(0);k<K;++k){
    draweta: etastart(k)  = (R::rnorm(eta(k),jumpeta(k))); 
    if(etastart(k)<0){goto draweta;}
    logCpetastart(k)      = logCpvMFcpp(p,etastart(k));
  }
  arma::mat logCpstart    = computelogCp(n, K, p, z, v, etastart, zeta);
  arma::mat loglikstart   = loglikelihood(n, K, logd, logb, logCpetastart, logCpzeta, logCpstart, Y);
  
  NumericVector temp1     = wrap((etastart - eta)/jumpeta);
  NumericVector temp2     = wrap((eta - etastart)/jumpeta);
  NumericVector lphi1     = Rcpp::pnorm(temp1, 0, 1, false, true);
  NumericVector lphi2     = Rcpp::pnorm(temp2, 0, 1, false, true);
  //Compute logalpha
  arma::rowvec logalpha   = (alphaeta-1.)*log(etastart/eta) + betaeta*(eta-etastart) + 
    (as<arma::rowvec>(lphi1) - as<arma::rowvec>(lphi2)) +  sum(loglikstart,0) - sum(loglik,0);
  arma::vec r_unif(K, arma::fill::randu);
  arma::uvec selupd       = arma::find(r_unif < exp(logalpha.t()));
  // update
  eta.elem(selupd)        = etastart.elem(selupd);
  logCp.cols(selupd)      = logCpstart.cols(selupd);
  loglik.cols(selupd)     = loglikstart.cols(selupd);
  logCpeta.elem(selupd)   = logCpetastart.elem(selupd);
  etaaccept.elem(selupd) += arma::ones<arma::rowvec>(selupd.n_elem); //Increase acceptance number to 1
}


//Update zeta
void zetaupdate(const double& n, const double& K, const double& p, const arma::mat& Y,
                const arma::mat& z, const arma::mat& v, const arma::vec& logd, const arma::rowvec& logb, const arma::rowvec& eta, double& zeta,
                const double& alphazeta, const double& betazeta, 
                const double& jumpzeta, double& zetaaccept, const arma::rowvec& logCpeta,
                double& logCpzeta, arma::mat& logCp, arma::mat& loglik){
  
  //Variable of the candidate value
  drawzeta: double zetastart=(R::rnorm(zeta,jumpzeta));
  
  //Draw from the proposal
  //if netative repeat
  if(zetastart<0){goto drawzeta;}
  
  //Compute Cpzetastart and phizetastart (see section 3.3)
  
  arma::mat logCpstart        = computelogCp(n, K, p, z, v, eta, zetastart);
  double logCpzetastart       = logCpvMFcpp(p,zetastart);
  arma::mat loglikstart       = loglikelihood(n, K, logd, logb, logCpeta, logCpzetastart, logCpstart, Y);
  
  double temp1          = (zetastart - zeta)/jumpzeta;
  double temp2          = (zeta - zetastart)/jumpzeta;
  double lphi1          = R::pnorm(temp1, 0, 1, false, true);
  double lphi2          = R::pnorm(temp2, 0, 1, false, true);
  double logalpha2zeta = (alphazeta-1.)*log(zetastart/zeta) + betazeta*(zeta-zetastart) + 
    lphi1 - lphi2 + accu(loglikstart) - accu(loglik);
  double logalphazeta=min(NumericVector::create(0,logalpha2zeta));
  
  if(unif_rand()<exp(logalphazeta)){
    zeta        = zetastart;      //Update zeta by zetastart
    logCp       = logCpstart;          //Update Cp by Cpstart
    loglik      = loglikstart;            //Update L by Lstart
    logCpzeta   = logCpzetastart;  //Update Cpzeta by Cpzetastart
    zetaaccept += 1;       //Increase acceptance number to 1
  }
}


///////////////////////////////////////////////////////////////////////////////
// The following function needs the initialization values and some 
// parameters of the model. It calls the previsous function to update   
// z, v, d, b, eta and zeta. it also updates the hypermarameters. Its arguments
// are described in section 3.4 of the technical details.
// [[Rcpp::export]]
List updateGP(const arma::mat& Y, const arma::mat& trait, const arma::mat& z0, const arma::mat& v0,
              const arma::vec& d0, const arma::rowvec& b0, const arma::rowvec& eta0, 
              const double& zeta0, const arma::uvec& fixv,
              const arma::uvec & consb, const double& nsimul, const bool& fdegrees, 
              const bool& fzeta, const NumericVector& hyperparms, const NumericVector& target,
              const NumericVector& jumpmin, const NumericVector& jumpmax, 
              const int& c, const bool& display_progress){
  
  //Parameters initionlization values
  arma::mat z = z0, v = v0;
  arma::vec d = d0;
  arma::vec logd = log(d);
  arma::rowvec b = b0, eta = eta0;
  arma::rowvec logb     = log(b);
  double zeta=zeta0 ;
  const double n=Y.n_rows;        // Nrow in ARD 
  const double K=Y.n_cols;        // Number of traits 
  const double p=z.n_cols;        // Hypersphere dimension 
  
  
  //Hyperparameters initionlization values
  
  double mud=hyperparms(0), sigmad=hyperparms(1);
  double mub=hyperparms(2), sigmab=hyperparms(3);
  double alphaeta=hyperparms(4), betaeta=hyperparms(5);
  double alphazeta=hyperparms(6), betazeta=hyperparms(7);
  
  //Jumping scales
  arma::vec jumpz = arma::ones(n), jumpd = arma::ones(n);
  arma::rowvec jumpb = arma::ones<arma::rowvec>(K), jumpeta = arma::ones<arma::rowvec>(K);
  double jumpzeta = 1;
  
  //Variables to count the number of accepance values
  arma::vec zaccept(n, arma::fill::zeros), daccept(n, arma::fill::zeros); 
  arma::rowvec baccept(K, arma::fill::zeros), etaaccept(K, arma::fill::zeros);
  double zetaaccept = 0;
  
  //Create output for parameters z, v, d, b, eta and zeta
  //We set the first entry equal to the initialization values
  List zoutput(nsimul+1); zoutput(0)=z;
  List voutput(nsimul+1); voutput(0)=v;
  arma::mat doutput(nsimul+1,n); doutput.row(0) = d.t();
  arma::mat boutput(nsimul+1,K); boutput.row(0) = b;
  arma::mat etaoutput(nsimul+1,K); etaoutput.row(0) = eta;
  NumericVector zetaoutput(nsimul+1); zetaoutput(0) = zeta;
  
  //Save hyperparameters updates
  //Note that alphaeta, betaeta, alphazeta and betazeta will not be updated
  NumericMatrix hyperupdate(nsimul+1,4);
  hyperupdate(0,_)=hyperparms;
  
  //Some useful variables
  //indexg is a list of size K. 
  //The entry k contains the indexg of individuals having the trait k
  //Cpeta initialization: (as in section 3.3). But here is the initization value
  //computed with the initialization value of eta.
  //Cpzeta initialization.
  List indexg(K);
  arma::rowvec logCpeta(K); 
  for(int k(0); k<K; ++k){
    logCpeta(k)=logCpvMFcpp(p,eta(k));
    indexg(k)=find(trait.col(k));
  }
  
  double logCpzeta=logCpvMFcpp(p,zeta);
  
  //Initialization of Cp
  arma::mat logCp  = computelogCp(n, K, p, z, v, eta, zeta);
  //Initialisation of L
  arma::mat loglik = loglikelihood(n, K, logd, logb, logCpeta, logCpzeta, logCp, Y);
  
  //save fixed values of v as vf (see section 3.3 update v). 
  const arma::mat vf=v.rows(fixv-1);
  //Compute PG0 (see section 3.3 update v)
  double maxlogb = max(logb.elem(consb-1));
  const double logPG0 = log(sum(exp(logb.elem(consb-1) - maxlogb))) + maxlogb;
  //Compute sums in Y rows and Y columns
  arma::colvec SumYrow=sum(Y,1);
  arma::rowvec SumYcol=sum(Y,0);
  
  //other variables
  double mubhat, sigmab2hat, mudhat, sigmad2hat, C0;
  
  //loop over nsimul
  Progress prgcpp(nsimul, display_progress);
  
  for(int t(1);t<nsimul+1;++t){
    //print step
    //if(t%100==1.0){std::cout<<"*";}
    //if(t%5000==0){std::cout<<" "<<floor(t/nsimul*100)<<"%"<<std::endl;}
    prgcpp.increment();
    
    //Step 2a update z
    zupdate(n, K, p, Y, z, v, logd, logb, eta, zeta, jumpz, zaccept, logCpeta, logCpzeta, logCp, loglik);
    
    //Step 2b update v
    vupdate(K, p, indexg, z, v, eta, fixv, vf);
    
    //Update Cp
    logCp = computelogCp(n, K, p, z, v, eta, zeta);
    
    //Step 2d update b
    bupdate(n, K, Y, z, v, logd, logb, eta, zeta, mub, sigmab, jumpb, baccept, logCpeta, logCpzeta, logCp, SumYcol);
    
    //Step 2c update d
    if(fdegrees==false){
      dupdate(n, K, Y, z, v, logd, logb, eta, zeta, mud, sigmad, jumpd, daccept, logCpeta, logCpzeta, logCp, SumYrow);
      
      //step 2e rescale d and b
      maxlogb = max(logb.elem(consb-1));
      C0      = log(sum(exp(logb.elem(consb-1) - maxlogb))) + maxlogb - logPG0;
      logd    = logd + C0;
      logb    = logb - C0;}
    
    //Update L
    loglik    = loglikelihood(n, K, logd, logb, logCpeta, logCpzeta, logCp, Y);
    
    //Step 2f update eta
    etaupdate(n, K, p, Y, z, v, logd, logb, eta, zeta, alphaeta, betaeta, jumpeta, etaaccept, logCpeta, logCpzeta, logCp, loglik);
    
    //Step 2g update zeta
    if(fzeta == false){
      zetaupdate(n, K, p, Y, z, v, logd, logb, eta, zeta, alphazeta, betazeta, jumpzeta, zetaaccept, logCpeta, logCpzeta, logCp, loglik);
    }
    
    //Step 2h update mub
    mubhat = mean(logb);
    mub    = (rnorm(1,mubhat,sigmab/sqrt(K)))(0);
    
    //step 2i update sigmab
    sigmab2hat = sum(pow(logb-mubhat,2));  
    sigmab     = pow(sigmab2hat/(rchisq(1,K-1)(0)),0.5);
    
    //Step 2j update mud
    mudhat = mean(logd);
    mud    = (rnorm(1,mudhat,sigmad/sqrt(n)))(0);
    
    //step 2k update sigmad
    sigmad2hat = sum(pow(logd-mudhat,2));  
    sigmad     = pow(sigmad2hat/(rchisq(1,n-1)(0)),0.5);
    
    //update jumping scale
    jumpz     -= (zaccept/t - target(0))/pow(t,c);
    fsetjump_v(jumpz, jumpmin(0), jumpmax(0));
    jumpd     += (daccept/t - target(1))/pow(t,c);
    fsetjump_v(jumpd, jumpmin(1), jumpmax(1));
    jumpb     += (baccept/t - target(2))/pow(t,c);
    fsetjump_r(jumpb, jumpmin(2), jumpmax(2));
    jumpeta   += (etaaccept/t - target(3))/pow(t,c);
    fsetjump_r(jumpeta, jumpmin(3), jumpmax(3));
    jumpzeta  += (zetaaccept/t - target(4))/pow(t,c);
    fsetjump_d(jumpzeta, jumpmin(4), jumpmax(4));
    
    //Save the updates
    zoutput(t)=z;
    voutput(t)= v;
    doutput.row(t)=exp(logd.t());
    boutput.row(t)=exp(logb);
    etaoutput.row(t)=eta;
    zetaoutput(t)=zeta;
    hyperupdate(t,0)=mud;
    hyperupdate(t,1)=sigmad;
    hyperupdate(t,2)=mub;
    hyperupdate(t,3)=sigmab;
  }
  
  //Create list for acceptance rate
  List acceptance=List::create(Named("z") = zaccept/nsimul, Named("d")   = daccept/nsimul,
                               Named("b") = baccept/nsimul, Named("eta") = etaaccept/nsimul,
                               Named("zeta") = zetaaccept/nsimul);
  // non updated hyper
  List hypernonupdate = List::create(Named("alphaeta")    = alphaeta,
                                      Named("betaeta")    = betaeta,
                                      Named("alphazeta")  = alphazeta,
                                      Named("betazeta")   = betazeta);
  
  // Save summary
  List simulations =  List::create(Named("z")    = zoutput,
                                   Named("v")    = voutput,
                                   Named("d")    = doutput,
                                   Named("b")    = boutput,
                                   Named("eta")  = etaoutput, 
                                   Named("zeta") = zetaoutput);
  
  
  List rethyperparms = List::create(Named("updated")     = hyperupdate,
                                    Named("non.updated") = hypernonupdate);
                                      
  return List::create(Named("simulations")       = simulations,
                      Named("hyperparms")        = rethyperparms,
                      Named("accept.rate")       = acceptance);
}

//////// Compute the graphs when Xnonard = NULL
// [[Rcpp::export]]
List dnetwork1(const double& T, const double& P, List& z, const arma::mat& d,
                    const arma::vec& zeta, const unsigned int& N,  const unsigned int& Metrostart,
                    const bool& display_progress){
  // Number of people with ARD
  const double ngraph=T-Metrostart+1;       // number of graphs
  double zetat, logCpzetat;
  arma::rowvec nuARDt(N);
  arma::rowvec nus(N, arma::fill::zeros), ds(N, arma::fill::zeros);
  
  
  arma::mat probt(N,N), prob(N,N,arma::fill::zeros), numat;
  
  //loop 
  Progress prgcpp(ngraph, display_progress);
  
  for(int t(0);t<ngraph;++t){  // thus convergence for t from Metrostart to T+1
    //print step
    //if((t+1)%100==1.0){std::cout<<"*";}
    //if((t+1)%5000==0){std::cout<<" "<<round((t+1)/ngraph*100)<<"%"<<std::endl;}
    prgcpp.increment();
    
    zetat           = zeta(t+Metrostart);          // extract zeta for itaration t+Metrostart
    arma::mat zt    = z(t+Metrostart);             // extract z for iteration t+Metrostart
    arma::rowvec dt = d.row(t+Metrostart);         // extract d for iteration t+Metrostart
    ds             += dt;
    logCpzetat      = logCpvMFcpp(P,zetat);        // Cp(P,zetat)
    //compute nu for ARD
    nuARDt          = log(dt) + 0.5*logCpzetat - 0.5*log(sum(dt)) ;
    nus            += nuARDt;
    
    // compute Probabilities
    numat = arma::repmat(nuARDt,N,1);
    probt = arma::exp(zetat*zt*zt.t() + numat + numat.t());
    probt.diag()=arma::zeros(N);   //zero on the diagonal
    probt*=(arma::sum(dt)/arma::accu(probt));
    
    prob += probt;
  }
  prob /= ngraph;
  ds   /= ngraph;
  nus  /= ngraph;
  
  arma::uvec tmp  = arma::find(prob > 1);
  prob.elem(tmp)  = arma::ones(tmp.n_elem);
  return  List::create(Named("dnetwork") = prob,
                       Named("degree")   = ds,
                       Named("nu")       = nus);
}



// [[Rcpp::export]]
List dnetwork2(const double& T, const double& P, List& z, const arma::mat& d,
                    const arma::vec& zeta, const arma::mat& Xard, const arma::mat& Xnonard, 
                    const arma::uvec& iARD, const arma::uvec& inonARD, const unsigned int& M, 
                    const unsigned int& Metrostart, const bool& display_progress){
  // Number of people with ARD
  const double n = Xard.n_rows;
  
  
  const double ngraph = T-Metrostart+1;       // number of graphs
  double zetat, logCpzetat;
  arma::rowvec nuARDt(n);
  
  // Number of people without ARD
  const double n2 = Xnonard.n_rows;
  //Number of peaple
  const double N  = n+n2;
  
  //compute neighbor and weight
  arma::mat neighbor(n2,M);            
  arma::mat weight(n2,M);
  cneighbor(n, n2, N, Xard, Xnonard, M, neighbor, weight);
  
  //Necessary variables
  arma::mat znonARDt(n2,P), ztall(N,P), probt(N,N), prob(N,N,arma::fill::zeros), numat;
  arma::uvec neighborj(M);
  arma::rowvec weightj(M), nunonARDt(n2), dnonARDt(n2), nut(N), nus(N, arma::fill::zeros),
  ds(N, arma::fill::zeros), dst(N);
  
  //loop 
  Progress prgcpp(ngraph, display_progress);
  
  for(int t(0);t<ngraph;++t){  // thus convergence for t from Metrostart to T+1
    //print step
    //if((t+1)%100==1.0){std::cout<<"*";}
    //if((t+1)%5000==0){std::cout<<" "<<round((t+1)/ngraph*100)<<"%"<<std::endl;}
    prgcpp.increment();
    
    zetat           = zeta(t+Metrostart);      // extract zeta for itaration t+Metrostart
    arma::mat zt    = z(t+Metrostart);        // extract z for iteration t+Metrostart
    arma::rowvec dt = d.row(t+Metrostart); // extract d for iteration t+Metrostart
    logCpzetat      = logCpvMFcpp(P,zetat);       // logCp(P,zetat)
    //compute nu for ARD
    nuARDt = log(dt) + 0.5*logCpzetat + 0.5*log(n/N) - 0.5*log(sum(dt)) ;
    
    //compute nu for non ARD
    for(int j(0);j<n2;++j){
      neighborj=arma::conv_to<arma::uvec>::from((neighbor.row(j)).t());
      weightj=weight.row(j);
      nunonARDt.col(j)=weightj*(nuARDt.elem(neighborj));
      znonARDt.row(j)=weightj*(zt.rows(neighborj));
    }
    znonARDt=normalise(znonARDt,2,1);
    
    // compute Probabilities
    ztall.rows(iARD)    = zt;
    ztall.rows(inonARD) = znonARDt;
    nut.elem(iARD)      = nuARDt;
    nut.elem(inonARD)   = nunonARDt;
    dst.elem(iARD)      = dt;

    numat = arma::repmat(nut,N,1);
    probt = arma::exp(zetat*ztall*ztall.t() + numat + numat.t());

    probt.diag()=arma::zeros(N); //zero on the diagonal
    
    // compute d for nonARD
    dnonARDt = (N/n)*arma::exp(nunonARDt)*sum(arma::exp(nuARDt))/exp(logCpzetat);
    probt   *= ((arma::sum(dt)+arma::sum(dnonARDt))/arma::accu(probt));
    
    dst.elem(inonARD) = dnonARDt;
    ds               += dst;
    nus              += nut;
    prob             += probt;
  }
  
  prob /= ngraph;
  ds   /= ngraph;
  nus  /= ngraph;
  
  arma::uvec tmp  = arma::find(prob > 1);
  prob.elem(tmp)  = arma::ones(tmp.n_elem);
  return  List::create(Named("dnetwork") = prob,
                       Named("degree")   = ds,
                       Named("nu")       = nus);
}