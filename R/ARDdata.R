#' @title Simulation of the distribution of the network
#' @description Compute the distribution of the network following  McCormick and Zheng (2015) and Breza et al. (2017).
#' @param nu is the vector of gregariousness.
#' @param d is the vector of degrees
#' @param zeta is a scale parameter that captures the influence of the latent positions on the link probabilities. 
#' @param z is a matrix where each row is a spherical coordinate.
#' @return a matrix of link probabilities
#' @references Breza, E., Chandrasekhar, A. G., McCormick, T. H., & Pan, M. (2017). Using aggregated relational data to feasibly
#'  identify network structure without network data (No. w23491). National Bureau of Economic Research. \url{https://arxiv.org/abs/1703.04157}.
#' @references McCormick, T. H., & Zheng, T. (2015). Latent surface models for networks using Aggregated Relational Data. 
#' Journal of the American Statistical Association, 110(512), 1684-1695. \url{https://www.tandfonline.com/doi/abs/10.1080/01621459.2014.991395}.
#' @examples 
#' \donotrun{
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
#' dist    <- sim.dnetwork(nu, d, zeta, z)}
#' @seealso \code{\link{sim.network}}
#' @export

sim.dnetwork  <- function(nu, d, zeta, z) {
  Prob(nu, d, zeta, z)
}

#' @title Simulation of the network
#' @description Generate networks from the network distribution.
#' @param dnetwork is the matrix of link probabilities.
#' @examples 
#' \donotrun{
#' N       <- 500 
#' zeta    <- 1
#' 
#' # Generate the spherical coordinates
#' z       <- rvMF(N, c(0, 0, 0))
#' 
#' # Genetate the gregariousness
#' nu      <- rnorm(N, -1.35, 0.37)
#' 
#' # Generale degrees
#' d       <- runif(N, 0, 45)
#' 
#' # Network distribution
#' dist    <- sim.dnetwork(nu, d, zeta, z)
#' 
#' # Generate network
#' G       <- sim.network(dist)}
#' @seealso \code{\link{sim.dnetwork}}
#' @export
sim.network <- function(dnetwork) {
  Graph(dnetwork)
}