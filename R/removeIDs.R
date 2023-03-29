#' @title Removes IDs with `NA` in a list of adjacency matrices in an optimal way.
#' @param network is a list of adjacency matrices
#' @param ncores is the number of cores to be used to run the program in parallel
#' @return List of adjacency matrices without missing values and a list of vectors of retained indeces
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%" 
#' @importFrom doRNG "%dorng%"
#' @examples 
#' A <- matrix(1:25, 5)
#' A[1, 1] <- NA
#' A[4, 2] <- NA
#' remove.ids(A)
#' 
#' B <- matrix(1:100, 10)
#' B[1, 1] <- NA
#' B[4, 2] <- NA
#' B[2, 4] <- NA
#' B[,8]   <-NA
#' remove.ids(B)
#' @export
remove.ids <- function(network, ncores = 1L){
  stopifnot(inherits(network, c("matrix", "data.frame", "list")))
  if(inherits(network, c("matrix", "data.frame"))) network <- list(network)
    
  # Construct cluster
  cl   <- makeCluster(ncores)
  registerDoParallel(cl)
  M    <- length(network)
  m    <- NULL
  out  <- foreach(m = 1:M, .packages  = "PartialNetwork") %dorng% {rem_non_fin(as.matrix(network[[m]]))}
  stopCluster(cl)
  network <- lapply(1:M, function(m) out[[m]]$net)
  id      <- lapply(1:M, function(m) c(out[[m]]$id))
  list(network = network, id = id)
}

