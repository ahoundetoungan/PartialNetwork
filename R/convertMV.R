#' @title Creating objects for network models
#' @description  `vec.to.mat` creates a list of square matrices from a given vector.
#' The elements of the generated matrices are taken from the vector and placed column-wise (ie. the first column is filled up before filling the second column) and from the first matrix of the list to the last matrix of the list. 
#' The diagonal of the generated matrices are zeros.
#' `mat.to.vec` creates a vector from a given list of square matrices .
#' The elements of the generated vector are taken from column-wise and from the first matrix of the list to the last matrix of the list, while dropping the diagonal entry.
#' `norm.network` row-normalizes matrices in a given list.
#' @param u numeric vector to convert.
#' @param W matrix or list of matrices to convert. 
#' @param N vector of sub-network sizes  such that `length(u) == sum(N*(N - 1))`.
#' @param normalise Boolean takes `TRUE` if the returned matrices should be row-normalized and `FALSE` otherwise.
#' @param ceiled Boolean takes `TRUE` if the given matrices should be ceiled before conversion and `FALSE` otherwise.
#' @param byrow Boolean takes `TRUE` is entries in the matrices should be taken by row and `FALSE` if they should be taken by column.
#' @return a vector of size `sum(N*(N - 1))` or list of `length(N)` square matrices. The sizes of the matrices are `N[1], N[2], ...`
#' @examples 
#' # Generate a list of adjacency matrices
#' ## sub-network size
#' N <- c(250, 370, 120)  
#' ## rate of friendship
#' p <- c(.2, .15, .18)   
#' ## network data
#' u <- unlist(lapply(1: 3, function(x) rbinom(N[x]*(N[x] - 1), 1, p[x])))
#' W <- vec.to.mat(u, N)
#' 
#' # Convert G into a list of row-normalized matrices
#' G <- norm.network(W)
#' 
#' # recover u
#' v <- mat.to.vec(G, ceiled = TRUE)
#' all.equal(u, v)
#' @seealso 
#' \code{\link{sim.network}}, \code{\link{sim.dnetwork}}, \code{\link{peer.avg}}.
#' @export
vec.to.mat <- function(u, N, normalise = FALSE, byrow = FALSE) {
  M        <- length(N)
  stopifnot(length(u) == sum(N*(N - 1)))
  out      <- NULL
  if (normalise) {
    out    <- frVtoMnorm(u, N, M)
  } else {
    out    <- frVtoM(u, N, M)
  }
  
  if(byrow) {
    out    <- lapply(out, t)
  }
  
  out
}


#' @rdname vec.to.mat
#' @export
mat.to.vec <- function(W, ceiled = FALSE, byrow = FALSE) {
  if (!is.list(W)) {
    if (is.matrix(W)) {
      W    <- list(W)
    } else {
      stop("W is neither a matrix nor a list")
    }
  }
  
  M        <- length(W)
  N        <- unlist(lapply(W, nrow))

  out      <- W
  if(byrow) {
    out    <- lapply(W, t)
  }
  if (ceiled) {
    out    <- frMceiltoV(out, N, M)
  } else {
    out    <- frMtoV(out, N, M)
  }
  
  out
}