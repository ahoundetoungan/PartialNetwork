#' @title Simulating Linear Peer Effect Models
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param parms A vector defining the true values of \eqn{(\lambda', \beta')'}, where  
#' \eqn{\lambda} is the peer effect parameter and \eqn{\beta} is the parameter vector of control variables (including contextual variables if any).
#' The parameters \eqn{\lambda} and \eqn{\beta} can also be specified separately using the arguments `lambda`, and `beta`.
#' @param lambda The true value of the vector \eqn{\lambda}.
#' @param beta The true value of the vector \eqn{\beta}.
#' @param epsilon A vector of idiosyncratic error terms. If not specified, it will be simulated from a standard normal distribution. 
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `simSAR` is called.
#' @param nthreads A strictly positive integer indicating the number of threads to use.
#' @description
#' `simSAR` simulates linear peer effect models.
#' @seealso \code{\link{smmSAR}}, \code{\link{mcmcSAR}}
#' @return A list containing:
#'     \item{y}{The simulated variable.}
#'     \item{Gy}{the average of y among friends.}
#' @examples 
#' set.seed(123)
#' ngr  <- 50
#' nvec <- rep(30, ngr)
#' n    <- sum(nvec)
#' G    <- lapply(1:ngr, function(z){
#'   Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z])
#'   diag(Gz) <- 0
#'   Gz/rowSums(Gz) # Row-normalized network
#' })
#' X    <- cbind(rnorm(n), rpois(n, 2))
#' l    <- 0.5
#' b    <- c(2, -0.5, 1)
#' 
#' out  <- simSAR(formula = ~ X, Glist = G, lambda = l, beta = b)
#' summary(out$y)
#' @export
#' @importFrom formula.tools env
#' @importFrom stats as.formula
#' @importFrom stats rnorm
simSAR   <- function(formula, Glist, parms, lambda, beta, epsilon, nthreads = 1, data){
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  igr      <- c(0, cumsum(nvec))
  # Data
  if (missing(data)) {
    data       <- env(formula)
  }
  formula      <- as.formula(formula)
  if(length(formula) != 2) stop("The `formula` argument is invalid. For data simulation, the expected format is `~ X1 + X2 + ...`.")
  
  ## call model.frame()
  mf          <- model.frame(formula, data = data)
  ## extract model matrices
  X           <- model.matrix(terms(formula, data = data, rhs = 1), mf)
  xname       <- colnames(X)
  
  Kx       <- ncol(X)
  eps      <- NULL
  if(missing(epsilon)){
    eps    <- rnorm(n)
  } else{
    eps    <- c(epsilon)
    if (!(length(eps) %in% c(1, n))) stop("`epsilon` must be either a scalar or an n-dimensional vector.")
    if (length(eps) == 1) eps <- rep(eps, n)
  }
  
  # parameters
  lamst    <- NULL
  lam      <- NULL
  b        <- NULL
  if (missing(parms)) {
    if (missing(lambda) | missing(beta)) {
      stop("Define either `parms` or `lambda` and `beta`.")
    }
    if (length(lambda) != 1){
      stop("lambda must be a scalar for the reduced-form model.")
    }
    lam   <- lambda
    
    if (length(beta) != Kx) stop("length(beta) is different from ncol(X).")
    b      <- beta
  } else{
    if (!missing(lambda) | !missing(beta)) {
      stop("Define either `parms` or `lambda` and `beta`.")
    }
    if (length(parms) != (1 + Kx)) stop("length(parms) is different from 1 + ncol(X).")
    lam   <- parms[1]
    b      <- tail(parms, Kx)
  }
  if (abs(lam) >= 1) {
    warning("The absolute value of peer effects is greater than or equal to one, which may lead to multiple or no equilibria.")
  }
  
  # Solving the game
  ## talpha
  talpha <- c(X %*% b + eps)
  y      <- fylim(G = Glist, talpha = talpha, igroup = igr, nvec = nvec, 
                  ngroup = M, lambda = lam, n = n, nthreads = nthreads)
  Gy     <- y$Gy
  y      <- y$y
  
  out    <- list("y"  = y,
                 "Gy" = Gy)
  out
}