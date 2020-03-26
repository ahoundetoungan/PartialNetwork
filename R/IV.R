#' @title Instrument Variables for SAR model
#' @description \code{sim.IV} generates Instrument Variables (IV) for linear-in-mean SAR models using only the distribution of the network. See Propositions 1 and 2 of Boucher and Houndetoungan (2020).
#' @param dnetwork is a square matrix (or a list of matrix if many groups) where the (i, j)th position is the probability of the event "i is linked to j".
#' @param X is a matrix of the individual observable characteristics.
#' @param y (optional) is the endogenous variable as a vector.
#' @param replication (optional, default = 1) is the number of repetitions (see details).
#' @param power (optional, default = 1) is the number of powers of the interaction matrix used to generate the instruments (see details).
#' 
#' @return A list of `replication` components. Each component is a list containing `G1y` (if the argument `y` was provided), `G1X` and `G2X` where `G1` and `G2` are independent draws of network from the distribution (see details).
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
#' @details Bramoullé et al. (2009) show that one can use \eqn{GX}, \eqn{G^2X}, ..., \eqn{G^P X} as instruments for \eqn{Gy}, where \eqn{P} is the maximal power desired.
#' \code{sim.IV} generate approximation of those instruments, based on Propositions 1 and 2 in Boucher and Houndetoungan (2020) (see also below).
#' The argument `power` is the maximal power desired.\cr
#' When \eqn{Gy} and the instruments \eqn{GX}, \eqn{G^2X}{GGX}, ..., \eqn{G^P X}{GG...GX} are not obseved, 
#' Boucher and Houndetoungan (2019) show that we can use one drawn from the distribution of the network in order to approximate \eqn{Gy}, but that
#' the same draw should not be used to approximate the instruments. Thus, each component in the function's output gives
#' `G1y` and `G1X` computed with the same network and `G2X` computed with another network, which can be used in order to approximate the instruments.
#' This process can be replicated several times and the argument `replication` can be used to set the number of replications desired.
#' @references 
#' Boucher, V., & Houndetoungan, A. (2020). Estimating peer effects using partial network data. \emph{Draft avaliable at} \url{https://houndetoungan.wixsite.com/aristide/research}.
#' @references 
#' Bramoullé, Y., Djebbari, H., & Fortin, B. (2009). Identification of peer effects through social networks. \emph{Journal of econometrics}, 150(1), 41-55. \url{https://www.sciencedirect.com/science/article/abs/pii/S0304407609000335}
#' @examples 
#' \dontrun{
#' library(AER)
#' # Number of groups
#' M             <- 5
#' # size of the group
#' N             <- 100      
#' # value of lambda (precision parameter for the network formation model)
#' lambda        <- 1   
#' # individual effects
#' beta          <- c(2, 1, 1.5) 
#' # contextual effects
#' gamma         <- c(5, -3) 
#' # endogenous effect
#' alpha         <- 0.4
#' # std-dev errors
#' se            <- 1
#' # heterogeneity of the linking probabilities
#' c             <- rnorm(N*N, 0, 1) 
#' 
#' # network probabilities
#' Prob          <- list()
#' 
#' for (m in 1:M) {
#'   # heterogeneity of the linking probabilities
#'   c           <- rnorm(N*N, 0, 1) 
#'   # generate linking probabilities
#'   Pm          <- matrix(exp(c / lambda) / (1 + exp(c / lambda)), N)
#'   # no self-link
#'   diag(Pm)    <- 0 
#'   Prob[[m]]   <- Pm
#' }
#' 
#' # generate data
#' X             <- matrix(data =  NA, nrow = 0, ncol = 2)
#' y             <- c()
#' 
#' for (m in 1:M) {
#'   # generate the 'observed network'
#'   Gm          <- sim.network(Prob[[m]]) 
#'   rs          <- rowSums(Gm)
#'   rs[rs == 0] <- 1
#'   # row-normalize
#'   Wm          <- Gm/rs
#'   # covariates
#'   Xm          <- cbind(rnorm(N,0,5),rpois(N,6)) 
#'   # endogenous variable, no contextual effect
#'   ym          <- solve(diag(N) - alpha * Wm) \%*\% (cbind(rep(1, N), Xm) \%*\% beta + rnorm(N,0,se)) 
#'   y           <- c(y, ym)
#'   X           <- rbind(X, Xm)
#' }
#' 
#' # generate instruments 
#' instr         <- sim.IV(Prob, X, y, replication = 1, power = 2)
#' 
#' GY1c1     <- instr[[1]]$G1y       # proxy for Gy (draw 1)
#' GXc1      <- instr[[1]]$G1X[,,1]  # proxy for GX (draw 1)
#' G2Xc1     <- instr[[1]]$G1X[,,2]  # proxy for GGX (draw 1)
#' GXc2      <- instr[[1]]$G2X[,,1]  # proxy for GX (draw 2)
#' G2Xc2     <- instr[[1]]$G2X[,,2]  # proxy for GGX (draw 2)
#' 
#' # build dataset
#' # keep only instruments constructed using a different draw than the one used to proxy Gy
#' dataset           <- as.data.frame(cbind(y,X,GY1c1,GXc2,G2Xc2)) 
#' # rename variables
#' colnames(dataset) <- c("y","X1","X2","Gy1","Z1","Z2","ZZ1","ZZ2")
#' results           <- ivreg(y ~ X1 + X2 + Gy1 | X1 + X2 + Z1 + Z2 + ZZ1 + ZZ2, data = dataset)
#' summary(results)
#' }
#' @seealso 
#' \code{\link{mcmcSAR}}
#' @importFrom abind abind
#' @export

sim.IV  <- function(dnetwork, X, y = NULL, replication = 1L, power = 1L) {
  
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
    out <- instruments2(dnetwork[[1]], X[r1:r2,], S = replication, pow = power)
  } else {
    out <- instruments1(dnetwork[[1]], X[r1:r2,], y[r1:r2], S = replication, pow = power)
  }
  
  if (M > 1)
  for (m in 2:M) {
    outm <- NULL
    r1    <- Ncum[m] + 1
    r2    <- Ncum[m + 1]
    if (is.null(y)) {
      outm <- instruments2(dnetwork[[m]], X[r1:r2,], S = replication, pow = power)
    } else {
      outm <- instruments1(dnetwork[[m]], X[r1:r2,], y[r1:r2], S = replication, pow = power)
    }
    
    for (s in 1:replication) {
      out[[s]]$G1y <- c(out[[s]]$G1y, outm[[s]]$G1y)
      out[[s]]$G1X <- abind(out[[s]]$G1X, outm[[s]]$G1X, along = 1)
      out[[s]]$G2X <- abind(out[[s]]$G2X, outm[[s]]$G2X, along = 1)
    }

  } 

  class(out) <- "sim.IV"
  out
}