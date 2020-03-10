#' @title Network Distribution
#' @description \code{fit.dnetwork} compute the network distribution from the simulations of the posterior distribution of the 
#' network formation model.
#' @param object is an `estim.ARD` object returned by \link{\code{mcmcARD}}
#' @param traitARD is the matrix of traits for individuals with ARD. The entry (i, k) is 1 if i has the trait k and 0 otherwise.
#' @param traitnonARD is the matrix of traits for individuals without ARD. The entry (j, k) is 1 if j has the trait k and 0 otherwise.
#' @param m is the number of neighbours used to compute the gregariousness and the degree for individuals without ARD.
#' @param burnin is the number of simulations from the posterior distribution used as burn-in. The network distribution will be computed
#' used the simulation from the iteration \code{burnin + 1}.
#' @param progress is logical; if TRUE, the progression will be printed in the console.
#' @details --- VINCENT Any details? ---
#' @return a matrix of the network distribution
#' @examples 
#' \donotrun{
#' FORTCOMING}
#' @export

fit.dnetwork   <- function(object, traitARD, traitnonARD = NULL, m, burnin = NULL, progress = TRUE){
  t1           <- Sys.time()
  if (class(object) != "estim.ARD") {
    stop("object class is not estim.ARD")
  }
  T          <- object$iteration
  P          <- object$p
  simu       <- object$simulations
  z          <- simu$z
  d          <- simu$d
  zeta       <- simu$zeta
  out        <- NULL
  Mst        <- round(T/2) + 1
  if (!is.null(burnin)){
    Mst      <- burnin + 1
  }
  
  if (Mst >= T) {
    stop("Burnin is too high")
  }
  
  if (is.null(traitnonARD)) {
    if (!missing(m)) {
      warning("m is defined but not used")
    }
    cat("ARD observed on the entire population \n")
    out      <- dnetwork1(T, P, z, d, zeta, traitARD, Mst, progress)
  } else {
    cat("ARD non observed on the entire population \n")
    out      <- dnetwork2(T, P, z, d, zeta, traitARD, traitnonARD, m, Mst, progress)
  }
  
  class(out) <- "dnetwork"
  
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs"))
  
  
  # Print the processing time
  cat("\n\n")
  cat("Average link probabilities estimated \n")
  cat("Iteration             :", T - Mst + 1, "\n  \n")
  
  nhours     <- floor(timer/3600)
  nminutes   <- floor((timer-3600*nhours)/60)%%60
  nseconds   <- timer-3600*nhours-60*nminutes
  cat("Elapsed time   : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n")
  
  out
}
