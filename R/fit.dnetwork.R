#' @title Fitting Network Distribution using ARD.
#' @description \code{fit.dnetwork} computes the network distribution using the simulations from the posterior distribution of the 
#' ARD network formation model. The linking probabilities are also computed for individuals without ARD if their traits are observed.
#' The degrees and the gregariousness of the individuals without ARD are computed from the sample with ARD using a m-nearest neighbors method.
#' @param object is an `estim.ARD` object returned by \code{\link{mcmcARD}}.
#' @param traitARD is the matrix of traits for individuals with ARD. The entry (i, k) is 1 if i has the trait k and 0 otherwise.
#' @param traitnonARD is the matrix of traits for individuals without ARD. The entry (j, k) is 1 if j has the trait k and 0 otherwise.
#' @param m is the number of neighbors used to compute the gregariousness and the degree for individuals without ARD.
#' @param burnin is the number of simulations from the posterior distribution used as burn-in. The network distribution will be computed
#' used the simulation from the iteration \code{burnin + 1}.
#' @param print is logical; if TRUE, the progression will be printed in the console.
#' @return a matrix of linking probabilities.
#' @examples 
#' \dontrun{
#' set.seed(123)
#' 
#' # GENERATE DATA
#' # Sample size
#' N  <- 500 
#' n  <- 300
#' 
#' # ARD parameters
#' genzeta <- 1
#' mu      <- -1.35
#' sigma   <- 0.37
#' K       <- 12    # number of traits
#' P       <- 3     # Sphere dimension 
#' 
#' # Generate z (spherical coordinates)
#' genz    <- rvMF(N,rep(0,P))
#' 
#' # Genetate nu  from a Normal distribution with parameters mu and sigma (The gregariousness)
#' gennu   <- rnorm(N,mu,sigma)
#' 
#' # compute degrees
#' gend <- N*exp(gennu)*exp(mu+0.5*sigma^2)*exp(logCpvMF(P,0) - logCpvMF(P,genzeta))
#' 
#' # Link probabilities
#' Probabilities <- sim.dnetwork(gennu,gend,genzeta,genz) 
#' 
#' # Adjacency matrix
#' G <- sim.network(Probabilities)
#' 
#' # Generate vk, the trait location
#' genv <- rvMF(K,rep(0,P))
#' 
#' # set fixed some vk  distant
#' genv[1,] <- c(1,0,0)
#' genv[2,] <- c(0,1,0)
#' genv[3,] <- c(0,0,1)
#' 
#' # eta, the intensity parameter
#' geneta   <-abs(rnorm(K,2,1))
#' 
#' # Build traits matrix
#' densityatz      <- matrix(0,N,K)
#' for(k in 1:K){
#'   densityatz[,k] <- dvMF(genz,genv[k,]*geneta[k])
#' }
#' 
#' trait       <- matrix(0,N,K)
#' for(k in 1:K){
#'   trait[,k] <- densityatz[,k]>sort(densityatz[,k],decreasing = T)[runif(1,0.05*N,0.25*N)]
#' }
#' 
#' # print a percentage of peaople having a trait
#' colSums(trait)*100/N
#' 
#' # Build ADR
#' ARD         <- G \%*\% trait
#' 
#' # generate b
#' genb        <- numeric(K)
#' for(k in 1:K){
#'   genb[k]   <- sum(G[,trait[,k]==1])/sum(G)
#' }
#' 
#' ############ ARD Posterior distribution ################### 
#' # EXAMPLE 1: ARD observed for the entire population
#' # initialization 
#' d0     <- exp(rnorm(N)); b0 <- exp(rnorm(K)); eta0 <- rep(1,K);
#' zeta0  <- 1; z0 <- matrix(rvMF(N,rep(0,P)),N); v0 <- matrix(rvMF(K,rep(0,P)),K)
#' 
#' # We need to fix some of the vk and bk for identification (see Breza et al. (2020) for details).
#' vfixcolumn      <- 1:5
#' bfixcolumn      <- c(3, 5)
#' b0[bfixcolumn]  <- genb[bfixcolumn]
#' v0[vfixcolumn,] <- genv[vfixcolumn,]
#' 
#' start  <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0, "zeta" = zeta0)
#' # MCMC ARD
#' out    <- mcmcARD(Y = ARD, traitARD = trait, start = start, fixv = vfixcolumn,
#'                   consb = bfixcolumn, iteration = 5000)
#'                   
#' # fit network distribution
#' dist   <- fit.dnetwork(out, traitARD = trait)
#' 
#' plot(rowSums(dist), gend)
#' 
#' # EXAMPLE 2: ARD observed for a sample of the population
#' # observed sample
#' traitard    <- trait[1:n, ]
#' traitnonard <- trait[(n + 1):N, ]
#' ARD         <- ARD[1:n,]
#' 
#' # initianalization 
#' d0     <- exp(rnorm(n)); b0 <- exp(rnorm(K)); eta0 <- rep(1,K);
#' zeta0  <- 1; z0 <- matrix(rvMF(n,rep(0,P)),n); v0 <- matrix(rvMF(K,rep(0,P)),K)
#' 
#' # We need to fix some of the vk and bk for identification (see Breza et al. (2020) for details).
#' vfixcolumn      <- 1:5
#' bfixcolumn      <- c(3, 5)
#' b0[bfixcolumn]  <- genb[bfixcolumn]
#' v0[vfixcolumn,] <- genv[vfixcolumn,]
#' 
#' start  <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0, "zeta" = zeta0)
#' # MCMC ARD
#' out    <- mcmcARD(Y = ARD, traitARD = traitard, start = start, fixv = vfixcolumn,
#'                   consb = bfixcolumn, iteration = 5000)
#'                   
#' # fit network distribution
#' dist   <- fit.dnetwork(out, traitARD = traitard, traitnonARD = traitnonard, m = 5)
#' 
#' plot(rowSums(dist)[1:n], gend[1:n], col = "blue")   # observed ard
#' points(rowSums(dist)[(n + 1):N], gend[(n + 1):N], col = "red")
#' }
#' @export

fit.dnetwork   <- function(object, traitARD, traitnonARD = NULL, m, burnin = NULL, print = TRUE){
  t1           <- Sys.time()
  if (class(object) != "estim.ARD") {
    stop("object class is not estim.ARD")
  }
  T          <- object$iteration
  P          <- object$p
  simu       <- object$simulations
  z          <- simu$z
  d          <- simu$d
  zeta       <- simu$zeta
  out        <- NULL
  Mst        <- round(T/2) + 1
  if (!is.null(burnin)){
    Mst      <- burnin + 1
  }
  
  if (Mst >= T) {
    stop("Burnin is too high")
  }
  
  if (is.null(traitnonARD)) {
    if (!missing(m)) {
      warning("m is defined but not used")
    }
    if (print) {
      cat("ARD observed on the entire population \n")
    }
    out      <- dnetwork1(T, P, z, d, zeta, traitARD, Mst, print)
  } else {
    if (print) {
      cat("ARD non observed on the entire population \n")
    }
    out      <- dnetwork2(T, P, z, d, zeta, traitARD, traitnonARD, m, Mst, print)
  }
  
  class(out) <- "dnetwork"
  
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs"))
  
  
  # Print the processing time
  if (print) {
    cat("\n\n")
    cat("Average link probabilities estimated \n")
    cat("Iteration             :", T - Mst + 1, "\n  \n")
    
    nhours     <- floor(timer/3600)
    nminutes   <- floor((timer-3600*nhours)/60)%%60
    nseconds   <- timer-3600*nhours-60*nminutes
    cat("Elapsed time   : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n")
  }
  out
}
