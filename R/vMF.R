#' @title Simulation from von Mises-Fisher Distribution
#' @description Random generation for the von Mises-Fisher Distribution
#' of dimension \code{p} with location parameter \code{mu} and intensity parameter \code{eta} 
#' (see Wood, 1994; Mardia, 2014).
#' @param size is the number of simulations.
#' @param theta is the parameter as \code{eta*mu}.
#' @importFrom Rcpp sourceCpp
#' @references 
#' Hornik, K., & Grün, B. (2014). \pkg{movMF}: An \R package for fitting mixtures of von Mises-Fisher distributions. \emph{Journal of Statistical Software}, 58(10), 1-31. \url{https://epub.wu.ac.at/4893/}.
#' @references  
#' Mardia, K. V. (2014). Statistics of directional data. Academic press. \url{https://www.elsevier.com/books/statistics-of-directional-data/mardia/978-0-12-471150-1}.
#' @references  
#' Wood, A. T. (1994). Simulation of the von Mises-Fisher distribution. Communications in statistics-simulation and computation, 23(1), 157-164. \url{https://www.tandfonline.com/doi/abs/10.1080/03610919408813161}.
#' @examples
#' # Draw 1000 vectors from vMF with parameters eta = 1 and mu = c(1,0)
#' rvMF(1000,c(1,0))
#' 
#' # Draw 10 vectors from vMF with parameters eta = sqrt(14) and mu proportional to (2,1,3)
#' rvMF(10,c(2,1,3))
#' 
#' # Draw from the vMF distribution with mean direction proportional to c(1, -1)
#' # and concentration parameter 3
#' rvMF(10, 3 * c(1, -1) / sqrt(2))
#' @export
rvMF  <- function(size, theta) {
  rvMFcpp(size, theta)
}



#' @title Normalization constant of von Mises-Fisher Distribution
#' @description log of the Normalization Constant for the von Mises-Fisher Distribution
#' of dimension \code{p} with location parameter \code{mu} and intensity parameter \code{eta}.
#' @param p is the dimension of the hypersphere.
#' @param eta is the intensity parameter.
#' @examples 
#' logCpvMF(2, 3.1)
#' @export
logCpvMF <- function(p, eta) {
  logCpvMFcpp(p, eta)
}

#' @title Density function of von Mises-Fisher Distribution
#' @description Density function for the von Mises-Fisher Distribution
#' of dimension \code{p} with location parameter equal to \code{mu} and intensity parameter \code{eta}.
#' @param z is a matrix where each row is a spherical coordinate at which the density will be evaluated.
#' @param theta is a vector of dimension `p` such that \code{eta*mu}.
#' @param log.p is logical; if TRUE, probabilities p are given as log(p).
#' @examples
#' # Draw 1000 vectors from vMF with parameter eta = 1 and mu = c(1,0)
#' z <- rvMF(1000, c(1,0))
#' 
#' # Compute the density at z
#' dvMF(z, c(1,0))
#' 
#' # Density of c(0, 1, 0, 0) with the parameter eta = 3 and mu = c(0, 1, 0, 0)
#' dvMF(matrix(c(0, 1, 0, 0), nrow = 1), c(0, 3, 0, 0)) 
#' @export
dvMF <- function(z, theta, log.p = FALSE) {
  dvMFcpp(z, theta, log.p)
}