#' @title Simulation of the distribution of the network for Breza et al. (2020)
#' @description Compute the distribution of the network following  McCormick and Zheng (2015) and Breza et al. (2020).
#' @param nu is the vector of gregariousness.
#' @param d is the vector of degrees.
#' @param zeta is a scale parameter that captures the influence of the latent positions on the link probabilities. 
#' @param z is a matrix where each row is a spherical coordinate.
#' @return a matrix of linking probabilities.
#' @references 
#' @references Breza, E., Chandrasekhar, A. G., McCormick, T. H., & Pan, M., 2020, Using aggregated relational data to feasibly identify network structure without network data, \emph{American Economic Review}, 110(8), 2454-84, \doi{10.1257/aer.20170861}
#' @references McCormick, T. H., & Zheng, T., 2015, Latent surface models for networks using Aggregated Relational Data, 
#' \emph{Journal of the American Statistical Association}, 110(512), 1684-1695, \doi{10.1080/01621459.2014.991395}.
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
