#' @title Bayesian Estimator of SAR model
#' 
#' 
#' @importFrom Formula as.Formula
#' @importFrom stats model.frame
#' @seealso \code{\link{sim.IV}}
#' @export
mcmcSAR <- function(formula,
                    contextual = TRUE,
                    start,
                    G0     = NULL,
                    hyperparms,
                    iteration = 2000,
                    ctrl.mcmc = list(),
                    data){
  t1           <- Sys.time()
  target       <- ctrl.mcmc$target
  jumpmin      <- ctrl.mcmc$jumpmin
  jumpmax      <- ctrl.mcmc$jumpmax
  print        <- ctrl.mcmc$print
  block.max    <- ctrl.mcmc$block.max
  c            <- ctrl.mcmc$c
  
  if (is.null(target)) {
    target     <- 0.44
  } else {
    if(length(target) != 1) {
      stop("target in ctrl.mcmc should be a scalar")
    }
  }
  
  if (is.null(jumpmin)) {
    jumpmin    <- 1e-12
  } else {
    if(length(jumpmin) != 1) {
      stop("jumpin in ctrl.mcmc should be a scalar")
    }
  }
  
  if (is.null(jumpmax)) {
    jumpmax    <- 10
  } else {
    if(length(jumpmax) != 1) {
      stop("jumpmax in ctrl.mcmc should be a scalar")
    }
  }
  
  if (is.null(block.max)) {
    block.max  <- 1
  } else {
    if(length(block.max) != 1) {
      stop("block.max in ctrl.mcmc should be a scalar")
    }
  }
  
  if (is.null(c)) {
    c          <- 0.6
  } else {
    if(length(c) != 1) {
      stop("c in ctrl.mcmc should be a scalar")
    }
  }
  
  f.t.data     <- formula.to.data(formula, contextual, data, ...)
  formula      <- f.t.data$formula
  Xone         <- f.t.data$Xone
  X            <- f.t.data$X
  y            <- f.t.data$y
  block.max    <- round(block.max)
  
  kbeta        <- ncol(Xone)
  kgamma       <- 0
  if (!is.null(X)) {
    kgamma     <- ncol(X)
  }
  
  col.xones   <- colnames(Xones)
  col.x       <- colnames(X)
  col.post    <- c(col.xones, paste0("G: ", col.x), "Gy", "se2")
  
  # hyperparameter
  dnetwork     <- hyperparms$dnetwork
  mutheta      <- hyperparms$mutheta
  invstheta    <- hyperparms$invstheta
  muzeta       <- hyperparms$muzeta
  invszeta     <- hyperparms$invszeta
  a            <- hyperparms$a
  b            <- hyperparms$b
  
  if (is.null(dnetwork)) {
    stop("dnetwork is not defined in hyperparms")
  }
  stopifnot(is.list(dnetwork))
  
  # M and network
  M            <- length(dnetwork)
  N            <- NULL
  if (is.null(X)) {
    if (is.null(G0)) {
      tmp      <- flistGnorm1nc(dnetwork, y, Xone, M)
      N        <- tmp$N
      G0       <- tmp$G
      y        <- tmp$ly
      Xone     <- tmp$lXone
    } else {
      tmp      <- flistGnorm2nc(dnetwork, G0, y, Xone, M)
      N        <- tmp$N
      G0       <- tmp$G
      y        <- tmp$ly
      Xone     <- tmp$lXone
    }
  } else {
    if (is.null(G0)) {
      tmp      <- flistGnorm1(dnetwork, y, Xone, X, M)
      N        <- tmp$N
      G0       <- tmp$G
      y        <- tmp$ly
      Xone     <- tmp$lXone
      X        <- tmp$lX
    } else {
      tmp      <- flistGnorm2(dnetwork, G0, y, Xone, X, M)
      N        <- tmp$N
      G0       <- tmp$G
      y        <- tmp$ly
      Xone     <- tmp$lXone
      X        <- tmp$lX
    }
  }
  
  
  
  if (is.null(X)) {
    stop("This version of the package does not allow MCMC for the model without contextual effect")
  } else {
    if (block.max == 1) {
      out      <- peerMCMC(y, X, Xone, G0, M, N, kbeta, kgamma, dnetwork, mutheta, invstheta, 
                            muzeta, invszeta, a, b, start, iteration, target, jumpmin, jumpmax,
                           c, print)
    } else {
      out      <- peerMCMCblock(y, X, Xone, G0, M, N, kbeta, kgamma, dnetwork, mutheta, invstheta, 
                                muzeta, invszeta, a, b, start, iteration, target, jumpmin, jumpmax,
                                c, block.max, print)
    }
  }
  
  
  colnames(out$posterior) <- col.post
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs")) 
  
  
  cat("\n\n")
  cat("The program successfully executed \n")
  cat("\n")
  cat("*************SUMMARY************* \n")
  cat("Number of group        : ", M, "\n")
  cat("Iteration              : ", iteration, "\n")
  
  # Print the processing time
  nhours     <- floor(timer/3600)
  nminutes   <- floor((timer-3600*nhours)/60)%%60
  nseconds   <- timer-3600*nhours-60*nminutes
  cat("Elapsed time           : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
  cat("Average acceptance rate: ", out$acceptance, "\n")

  
  return(c(list("n.group" = M, "N" = c(N), "iteration" = iteration),
                out, 
                list("hyperparms" = hyperparms, "formula" = formula)))
}