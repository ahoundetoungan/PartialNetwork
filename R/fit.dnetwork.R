#' @title Fitting Network Distribution using ARD.
#' @description \code{fit.dnetwork} computes the network distribution using the simulations from the posterior distribution of the 
#' ARD network formation model. The linking probabilities are also computed for individuals without ARD.
#' The degrees and the gregariousness of the individuals without ARD are computed from the sample with ARD using a k-nearest neighbors method.
#' @param object `estim.ARD` object returned by \code{\link{mcmcARD}}.
#' @param X (required when ARD are available for a sample of individuals) is a matrix of variables describing individuals with ARD and those without ARD. This matrix will be used to compute distance between individuals
#' in the k-nearest neighbors approach.
#' This could be the matrix of traits (see details).
#' @param obsARD logical vector of length `nrow(X)` (number of individuals with and without ARD), where the i-th entry equal to `TRUE` if the i-th individual in `X` has ARD and `FALSE` otherwise.
#' If missing, `obsARD = rep(c(TRUE, FALSE), n1, n2)`, where `n1` is the number of individuals with ARD (see details).
#' @param m number of neighbors used to compute the gregariousness and the degree for individuals without ARD (default value is `1`).
#' @param burnin number of simulations from the posterior distribution used as burn-in. The network distribution will be computed
#' used the simulation from the iteration \code{burnin + 1}.
#' @param print logical; if TRUE, the progression will be printed in the console.
#' @return A list consisting of:
#'         \item{dnetwork}{posterior mean of the network distribution.} 
#'         \item{degree}{posterior mean of the degree.} 
#'         \item{nu}{posterior mean of the gregariousness, nu.} 
#' @details The order of individuals provided through the arguments `traitARD` and `ARD` (when calling the function \code{\link{mcmcARD}}) should fit the order of individuals in
#' `X` and `obsARD`. Especially, the i-th row of `X[obsARD,]` should correspond to the i-th row in `traitARD` or `ARD`.
#' @examples 
#' set.seed(123)
#' # GENERATE DATA
#' # Sample size
#' N  <- 50
#' n  <- 30
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
#' NK          <- floor(runif(K, 0.8, 0.95)*colSums(densityatz)/apply(densityatz, 2, max))
#' for (k in 1:K) {
#'   trait[,k]       <- rbinom(N, 1, NK[k]*densityatz[,k]/sum(densityatz[,k]))
#' }
#' 
#' # print a percentage of people having a trait
#' colSums(trait)*100/N
#' 
#' # Build ARD
#' ARD         <- G %*% trait
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
#' vfixcolumn      <- 1:6
#' bfixcolumn      <- c(3, 5)
#' b0[bfixcolumn]  <- genb[bfixcolumn]
#' v0[vfixcolumn,] <- genv[vfixcolumn,]
#' 
#' start  <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0, "zeta" = zeta0)
#' # # MCMC ARD
#' # REMOVE COMMENTS TO RUN THE REST OF THE CODE.
#' # mcmcARD USES Rcpp FUNCTIONS IN PARALLEL TO BE FAST AND CRAN POLICY DOES 
#' # NOT ALLOW CODE IN PARALLEL. FOR THIS REASON WE PUT THE REST OF THE CODE IN 
#' # PARALLEL SO THAT CRAN ACCEPTS THE PACKAGE.
#' # out    <- mcmcARD(Y = ARD, traitARD = trait, start = start, fixv = vfixcolumn,
#' #                   consb = bfixcolumn, iteration = 500)
#' # 
#' # # fit network distribution
#' # dist   <- fit.dnetwork(out)
#' # 
#' # plot(rowSums(dist$dnetwork), gend)
#' # abline(0, 1, col = "red")
#' # 
#' # # EXAMPLE 2: ARD observed for a sample of the population
#' # # observed sample
#' # selectARD   <- sort(sample(1:N, n, FALSE))
#' # traitard    <- trait[selectARD,]
#' # ARD         <- ARD[selectARD,]
#' # logicalARD  <- (1:N) %in% selectARD
#' # 
#' # # initianalization
#' # d0     <- exp(rnorm(n)); b0 <- exp(rnorm(K)); eta0 <- rep(1,K);
#' # zeta0  <- 1; z0 <- matrix(rvMF(n,rep(0,P)),n); v0 <- matrix(rvMF(K,rep(0,P)),K)
#' # 
#' # # We need to fix some of the vk and bk for identification (see Breza et al. (2020) for details).
#' # vfixcolumn      <- 1:6
#' # bfixcolumn      <- c(3, 5)
#' # b0[bfixcolumn]  <- genb[bfixcolumn]
#' # v0[vfixcolumn,] <- genv[vfixcolumn,]
#' # 
#' # start  <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0, "zeta" = zeta0)
#' # # MCMC ARD
#' # REMOVE COMMENTS TO RUN THE REST OF THE CODE.
#' # mcmcARD USES Rcpp FUNCTIONS IN PARALLEL TO BE FAST AND CRAN POLICY DOES 
#' # NOT ALLOW CODE IN PARALLEL. FOR THIS REASON WE PUT THE REST OF THE CODE IN 
#' # PARALLEL SO THAT CRAN ACCEPTS THE PACKAGE.
#' # out    <- mcmcARD(Y = ARD, traitARD = traitard, start = start, fixv = vfixcolumn,
#' #                   consb = bfixcolumn, iteration = 500)
#' # 
#' # # fit network distribution
#' # dist   <- fit.dnetwork(out, X = trait, obsARD = logicalARD, m = 1)
#' # 
#' # library(ggplot2)
#' # ggplot(data.frame("etimated.degree" = dist$degree,
#' #                   "true.degree"     = gend,
#' #                   "observed"        = ifelse(logicalARD, TRUE, FALSE)),
#' #        aes(x = etimated.degree, y = true.degree, colour = observed)) +
#' #   geom_point()
#' @export

fit.dnetwork   <- function(object, X = NULL, obsARD = NULL,
                           m = NULL, burnin = NULL, print = TRUE){
  t1         <- Sys.time()
  stopifnot(inherits(object, "estim.ARD"))
  T          <- object$iteration
  P          <- object$p
  simu       <- object$simulations
  N1         <- object$n
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
  
  if (is.null(X)) {
    if (print) {
      cat("ARD observed on the entire population \n")
    }
    out      <- dnetwork1(T, P, z, d, zeta, N1, Mst, print)
  } else {
    if (print) {
      cat("ARD non observed on the entire population \n")
    }
    stopifnot(is.matrix(X))
    N        <- nrow(X)
    N2       <- N - N1
    if(N <= object$n) {
      stop("X provided when nrow(X) is not greater than the ARD sample size")
    }
    if(!is.null(obsARD)) {
      stopifnot(is.vector(obsARD))
      stopifnot(is.logical(obsARD))
      stopifnot(length(obsARD) == nrow(X))
    } else {
      obsARD <- rep(c(TRUE, FALSE), c(N1, N2))
    }
    iARD     <- which(obsARD)
    inonARD  <- which(!obsARD)
    XARD     <- X[iARD,]
    XnonARD  <- X[inonARD,]
    iARD     <- iARD - 1
    inonARD  <- inonARD - 1
    if(is.null(m)) {
      m      <- 1
    }
      
    out      <- dnetwork2(T, P, z, d, zeta, XARD, XnonARD, iARD, inonARD, m, Mst, print)
  }
  out$degree <- c(out$degree)
  out$nu     <- c(out$nu)
  
  #class(out) <- "dnetwork"
  
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs"))
  
  
  # Print the processing time
  if (print) {
    cat("\n")
    cat("Average link probabilities estimated \n")
    cat("Iteration      :", T - Mst + 1, "\n")
    
    nhours     <- floor(timer/3600)
    nminutes   <- floor((timer-3600*nhours)/60)%%60
    nseconds   <- timer-3600*nhours-60*nminutes
    cat("Elapsed time   : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n")
  }
  out
}
