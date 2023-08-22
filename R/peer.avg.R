#' @title Computing peer average value
#' @description `peer.avg` computes peer average value using network data (as a list) and observable characteristics.
#' @param Glist the adjacency matrix or list sub-adjacency matrix.
#' @param V vector or matrix of observable characteristics.
#' @param export.as.list (optional) boolean to indicate if the output should be a list of matrices or a single matrix.
#' @return the matrix product `diag(Glist[[1]], Glist[[2]], ...) %*% V`, where `diag()` is the block diagonal operator.
#' @examples 
#' # Generate a list of adjacency matrices
#' ## sub-network size
#' N  <- c(250, 370, 120)  
#' ## rate of friendship
#' p  <- c(.2, .15, .18)   
#' ## network data
#' u  <- unlist(lapply(1: 3, function(x) rbinom(N[x]*(N[x] - 1), 1, p[x])))
#' G  <- vec.to.mat(u, N, normalise = TRUE)
#' 
#' # Generate a vector y
#' y  <- rnorm(sum(N))
#' 
#' # Compute G%*%y
#' Gy <- peer.avg(Glist = G, V = y)
#' @seealso 
#' \code{\link{sim.network}}
#' @export



peer.avg    <- function(Glist, V, export.as.list = FALSE) {
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  
  v.is.mat  <- !is.null(dim(V))
  V         <- as.matrix(V)
  
  cnames    <- colnames(V)
  if(!is.null(cnames)) {
    cnames  <-  paste0("G.", cnames)
  }
  
  M         <- length(Glist)
  N         <- unlist(lapply(Glist, ncol))
  
  if (sum(N) != nrow(V)) {
    stop("Glist and V do not match")
  }
  Ncum      <- c(0, cumsum(N))
  
  out       <- lapply(1:M, peer.avg.single, Glist = Glist, V = V, Ncum = Ncum, cnames = cnames, v.is.mat = v.is.mat)
  if (!export.as.list) {
    out     <- do.call(rbind, out)
  }
  if(!v.is.mat) out <- c(out)
  out
}


peer.avg.single <- function(m, Glist, V, Ncum, cnames, v.is.mat) {
  if (!is.matrix(Glist[[m]])) {
    stop("All components in Glist are not matrices")
  }
  out <- Glist[[m]]%*%V[(Ncum[m] + 1):Ncum[m + 1],,drop = FALSE]
}
