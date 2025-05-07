#' @title Instrument Variables for SAR model
#' @description \code{sim.IV} generates Instrument Variables (IV) for linear-in-mean SAR models using only the distribution of the network. See Propositions 1 and 2 of Boucher and Houndetoungan (2020).
#' @param dnetwork network matrix of list of sub-network matrices, where the (i, j)-th position is the probability that i be connected to j.
#' @param X matrix of the individual observable characteristics.
#' @param y (optional) the endogenous variable as a vector.
#' @param replication (optional, default = 1) is the number of repetitions (see details).
#' @param power (optional, default = 1) is the number of powers of the interaction matrix used to generate the instruments (see details).
#' @param exp.network (optional, default = FALSE) indicates if simulated network should be exported.
#' 
#' @return list of `replication` components. Each component is a list containing `G1y` (if the argument `y` was provided), `G1` (if `exp.network = TRUE`), `G2` (if `exp.network = TRUE`) , `G1X`, and 
#' `G2X` where `G1` and `G2` are independent draws of network from the distribution (see details).
#' \item{G1y}{is an approximation of \eqn{Gy}.}
#' \item{G1X}{is an approximation of \eqn{G^pX}{GG ... GX}
#' with the same network draw as that used in `G1y`. `G1X` is an array of dimension \eqn{N \times K \times power}{N * K * power}, where \eqn{K}{K} is the number of column in 
#' `X`. For any \eqn{p \in \{1, 2, ..., power\}}{p = 1, 2, ..., power}, the approximation of \eqn{G^pX}{GG ... GX}
#' is given by \code{G1X[,,p]}.}
#' \item{G2X}{is an approximation of \eqn{G^pX}{GG ... GX}
#' with a different different network. `G2X` is an array of dimension \eqn{N \times K \times power}{N * K * power}. 
#' For any \eqn{p \in \{1, 2, ..., power\}}{p = 1, 2, ..., power}, the approximation of \eqn{G^pX}{GG ... GX}
#' is given by \code{G2X[,,p]}.}
#' 
#' 
#' @details Bramoulle et al. (2009) show that one can use \eqn{GX}, \eqn{G^2X}, ..., \eqn{G^P X} as instruments for \eqn{Gy}, where \eqn{P} is the maximal power desired.
#' \code{sim.IV} generate approximation of those instruments, based on Propositions 1 and 2 in Boucher and Houndetoungan (2020) (see also below).
#' The argument `power` is the maximal power desired.\cr
#' When \eqn{Gy} and the instruments \eqn{GX}, \eqn{G^2X}{GGX}, ..., \eqn{G^P X}{GG...GX} are not observed, 
#' Boucher and Houndetoungan (2022) show that we can use one drawn from the distribution of the network in order to approximate \eqn{Gy}, but that
#' the same draw should not be used to approximate the instruments. Thus, each component in the function's output gives
#' `G1y` and `G1X` computed with the same network and `G2X` computed with another network, which can be used in order to approximate the instruments.
#' This process can be replicated several times and the argument `replication` can be used to set the number of replications desired.
#' @examples 
#' library(AER)
#' # Number of groups
#' M             <- 30
#' # size of each group
#' N             <- rep(50,M)
#' # individual effects
#' beta          <- c(2,1,1.5) 
#' # endogenous effects
#' alpha         <- 0.4
#' # std-dev errors
#' se            <- 2 
#' # prior distribution
#' prior         <- runif(sum(N*(N-1)))
#' prior         <- vec.to.mat(prior, N, normalise = FALSE)
#' # covariates
#' X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
#' # true network
#' G0            <- sim.network(prior)
#' # normalise 
#' G0norm        <- norm.network(G0)
#' # simulate dependent variable use an external package
#' y             <- CDatanet::simsar(~ X, Glist = G0norm,
#'                                   theta = c(alpha, beta, se))
#' y             <- y$y
#' # generate instruments 
#' instr         <- sim.IV(prior, X, y, replication = 1, power = 1)
#' 
#' GY1c1         <- instr[[1]]$G1y       # proxy for Gy (draw 1)
#' GXc1          <- instr[[1]]$G1X[,,1]  # proxy for GX (draw 1)
#' GXc2          <- instr[[1]]$G2X[,,1]  # proxy for GX (draw 2)
#' # build dataset
#' # keep only instrument constructed using a different draw than the one used to proxy Gy
#' dataset           <- as.data.frame(cbind(y, X, GY1c1, GXc1, GXc2)) 
#' colnames(dataset) <- c("y","X1","X2","G1y", "G1X1", "G1X2", "G2X1", "G2X2") 
#' 
#' # Same draws
#' out.iv1           <- ivreg(y ~ X1 + X2 + G1y | X1 + X2 + G1X1 + G1X2, data = dataset)
#' summary(out.iv1)
#' 
#' # Different draws
#' out.iv2           <- ivreg(y ~ X1 + X2 + G1y | X1 + X2 + G2X1 + G2X2, data = dataset)
#' summary(out.iv2)
#' @seealso 
#' \code{\link{mcmcSAR}}
#' @importFrom abind abind
#' @export

sim.IV  <- function(dnetwork, X, y = NULL, replication = 1L, power = 1L, exp.network = FALSE) {
  
  M     <- NULL
  if (is.list(dnetwork)) {
    M          <- length(dnetwork)
  } else {
    if (is.matrix(dnetwork)) {
      M        <- 1
      dnetwork <- list(dnetwork)
    } else {
      stop("dnetwork is neither a matrix nor a list")
    }
  }
  
  if (! is.matrix(X)) {
    X   <- as.matrix(X)
  }
  out   <- NULL
  Nvec  <- unlist(lapply(dnetwork, nrow))
  Ncum  <- c(0, cumsum(Nvec))
  r1    <- Ncum[1] + 1
  r2    <- Ncum[2]
  if (is.null(y)) {
    out <- instruments2(dnetwork[[1]], X[r1:r2,], S = replication, pow = power, expG = exp.network)
  } else {
    out <- instruments1(dnetwork[[1]], X[r1:r2,], y[r1:r2], S = replication, pow = power, expG = exp.network)
  }
  
  if (M > 1)
  for (m in 2:M) {
    outm <- NULL
    r1    <- Ncum[m] + 1
    r2    <- Ncum[m + 1]
    if (is.null(y)) {
      outm <- instruments2(dnetwork[[m]], X[r1:r2,], S = replication, pow = power, expG = exp.network)
    } else {
      outm <- instruments1(dnetwork[[m]], X[r1:r2,], y[r1:r2], S = replication, pow = power, expG = exp.network)
    }
    
    for (s in 1:replication) {
      out[[s]]$G1y  <- c(out[[s]]$G1y, outm[[s]]$G1y)
      out[[s]]$G1X  <- abind(out[[s]]$G1X, outm[[s]]$G1X, along = 1)
      out[[s]]$G2X  <- abind(out[[s]]$G2X, outm[[s]]$G2X, along = 1)
      if(exp.network) {
        out[[s]]$G1 <- c(out[[s]]$G1, outm[[s]]$G1)
        out[[s]]$G2 <- c(out[[s]]$G2, outm[[s]]$G2)
      }
    }

  } 

  class(out) <- "sim.IV"
  out
}