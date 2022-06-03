#' @title Simulation of the distribution of the network for Breza et al. (2020)
#' @description Compute the distribution of the network following  McCormick and Zheng (2015) and Breza et al. (2020).
#' @param nu is the vector of gregariousness.
#' @param d is the vector of degrees.
#' @param zeta is a scale parameter that captures the influence of the latent positions on the link probabilities. 
#' @param z is a matrix where each row is a spherical coordinate.
#' @return a matrix of linking probabilities.
#' @examples 
#' N       <- 500 
#' zeta    <- 1
#' 
#' # Generate the spherical coordinates
#' z       <- rvMF(N, c(0, 0, 0))
#' 
#' # Genetate the gregariousness
#' nu      <- rnorm(N, -1.35, 0.37)
#' 
#' # Generate degrees
#' d       <- runif(N, 0, 45)
#' 
#' dist    <- sim.dnetwork(nu, d, zeta, z)
#' @seealso \code{\link{sim.network}}
#' @export

sim.dnetwork  <- function(nu, d, zeta, z) {
  Prob(nu, d, zeta, z)
}
