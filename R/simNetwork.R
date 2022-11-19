#' @title Simulating network data
#' @param dnetwork is a list of sub-network matrices, where the (i, j)-th position of the m-th matrix is the probability that i be connected to j, with i and j individuals from the m-th network.
#' @param normalise boolean takes `TRUE` if the returned matrices should be row-normalized and `FALSE` otherwise.
#' @return list of (row-normalized) adjacency matrices.
#' @examples 
#' # Generate a list of adjacency matrices
#' ## sub-network size
#' N         <- c(250, 370, 120)  
#' ## distribution
#' dnetwork  <- lapply(N, function(x) matrix(runif(x^2), x))
#' ## network
#' G         <- sim.network(dnetwork)
#' @seealso \code{\link{sim.dnetwork}}
#' @export
sim.network <- function(dnetwork, normalise = FALSE) {
  trsf          <- FALSE
  if (!is.list(dnetwork)) {
    if (is.matrix(dnetwork)) {
      dnetwork  <- list(dnetwork)
      trsf      <- TRUE
    } else {
      stop("dnetwork is neither a matrix nor a list of matrices")
    }
  }

  M        <- length(dnetwork)
  N        <- unlist(lapply(dnetwork, nrow))
  out      <- NULL
  if (normalise) {
    out    <- simGnorm(dnetwork = dnetwork, N = N, M = M)
  } else {
    out    <- simG(dnetwork = dnetwork, N = N, M = M)
  }
  if (trsf) {
    out    <- out[[1]]
  }
  out
}


#' @rdname vec.to.mat
#' @export
norm.network <- function(W) {
  trsf     <- FALSE
  if (!is.list(W)) {
    if (is.matrix(W)) {
      W    <- list(W)
      trsf <- TRUE
    } else {
      stop("W is neither a matrix nor a list of matrices")
    }
  }
  stopifnot(inherits(W, "list"))

  M        <- length(W)
  out      <- fGnormalise(W, M)
  if (trsf) {
    out    <- out[[1]]
  }
  out
}

