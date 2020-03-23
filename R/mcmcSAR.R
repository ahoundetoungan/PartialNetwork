#' @title Bayesian Estimator of SAR model
#' @description \code{mcmcSAR} implements the Bayesian estimator of the linear-in-mean SAR model when only the linking probabilities are available.
#' The linking probabilities are used as hyperparameter of the prior distribution of the network.
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' the `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param  start (optional) is the vector of starting value of the model parameter as \eqn{(\beta' ~ \gamma' ~ \alpha ~ \sigma^2)'}{(\beta'  \gamma'  \alpha  se2)'},
#' where \eqn{\beta} is the individual variables parameter, \eqn{\gamma} is the contextual variables parameter, \eqn{\alpha} is the peer effect parameter
#' and \eqn{\sigma^2}{se2} the variance of the error term. If the `start` is missing, a Maximum Likelihood estimator will be used, where
#' the newotk matrix is that given through the argument `G0` (if provided) or generated from it distribution `dnetwork` (see argument `hyperparms`).
#' @param hyperparms is a list of hyperparameter parameter. It should contain at least `dnetwork`, the linking probabilities. As there are `M` groups and
#' individual from different groups are not linked, `dnetwork` is a list of `M` matrix where each matrix of the link probabilities in one group.
#' @param G0 (optional) is the starting value of the (row normalized) network as a list of `M` sub-network.
#' @param iteration is the number of MCMC steps to be performed.
#' @param ctrl.mcmc is a list of MCMC controls (See details).
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @details The model is given by
#' \deqn{\mathbf{y} = \mathbf{X}\beta + \mathbf{G}\mathbf{X}\gamma + \alpha \mathbf{G}\mathbf{y} + \epsilon.}{y = X\beta + GX\gamma + \alpha Gy + \epsilon.}
#' The parameters to estimate in this model are the matrix \eqn{\mathbf{G}}{G}, the vectors \eqn{\beta}, \eqn{\gamma} and the scalar \eqn{\alpha}, \eqn{\sigma^2}{se2}.
#' Prior distributions are assumed on \eqn{\mathbf{A}}, the adjacency matrix in which \eqn{\mathbf{A}_{ij} = 1}{A[i,j] = 1} if i is  connected to j and
#' \eqn{\mathbf{A}_{ij} = 0}{A[i,j] = 0} otherwise, and on \eqn{\beta}, \eqn{\gamma}, \eqn{\alpha} and \eqn{\sigma^2}{se2}.
#' \deqn{\mathbf{A}_{ij} \sim Bernoulli(\mathbf{P}_{ij})}{A[i,j] ~ Bernoulli(P[i,j])}
#' \deqn{(\beta' ~ \gamma')'|\sigma^2 \sim \mathcal{N}(\mu_{\theta}, \sigma^2\Sigma_{\theta})}{(\beta' ~ \gamma')'|se2 ~ N(mutheta, se2*stheta)}
#' \deqn{\zeta = \log\left(\frac{\alpha}{1 - \alpha}\right) \sim \mathcal{N}(\mu_{\zeta}, \sigma_{\zeta})}{\zeta = log(\alpha/(1 - \alpha)) ~ N(muzeta, szeta)}
#' \deqn{\sigma^2 \sim IG(\frac{a}{2}, \frac{b}{2})}{se2 ~ IG(a/2, b/2)}
#' where \eqn{\mathbf{P}}{P} is the linking probability.\cr
#' 
#' All the hyperparametera can be defined through the argument `hyperparms` (a list)  and should be named as follow.
#' \itemize{
#' \item `dnetwork`, the linking probabilities (list of M matrix). This hyperparameter is required.
#' \item `mutheta`, the prior mean of \eqn{(\beta' ~ \gamma')'|\sigma^2}{(\beta' ~ \gamma')'|se2}. The default value assumes that
#' the prior mean is zero.
#' \item `invstheta` as \eqn{\Sigma_{\theta}^{-1}}{inverse of `stheta`}. The default value is a diagonal matrix with 0.01 on the diagonal.
#' \item `muzeta`, the prior mean of \eqn{\zeta}. The default value is zero.
#' \item `invszeta`, the inverse of the prior variance of \eqn{\zeta} with default value equal to 2.
#' \item `a` and `b` which default values equal to 4.2 and 2.2 respectively. This means for example that the prior mean of \eqn{\sigma^2}{se2} is 1.
#' }
#' Inverses are used for the priori variance through the argument `hyperparms`  in order to allow non informative prior. Set the inverse of the prior
#' variance to 0 is equivalent to assume a non informative prior.\cr
#' During the MCMC, the jumping scale of \eqn{\alpha} is updated following Atchadé and Rosenthal (2005) in order to target the acceptance rate to the `target` value. This
#' requires to set a minimal and a maximal jumping scales through the parameter `ctrl.mcmc`. The parameter `ctrl.mcmc` is a list which can contain the following named components.
#' \itemize{
#' \item{`target`}: the default value is 0.44. 
#' \item{`jumpmin`}: the default value is \code{1e-12}. 
#' \item{`jumpmax`}: the default value is \code{10}. 
#' \item{`print.level`}: an integer in \{0, 1, 2\} that indicates if the MCMC progression should be printed in the console.
#'  If 0, the MCMC progression is not be printed. If 1 (default value), the progression is printed and if 2,
#'  the simulations from the posterior distribution are printed.
#' \item{`block.max`}: The maximal number of entries that can be updated simultaneously in \eqn{\mathbf{A}}{A}. For some convergence reasons, it mighit be 
#' more efficient to update simultaneously 2 or 3 entries (see Boucher and Houndetoungan, 2019). 
#' }
#' If `block.max` > 1, several entries are ramdomly chosen from the same row and updated simultaneously. The number of entries chosen is randomly 
#' chosen between 1 and `block.max`. In addition, the entries are not chosen in order. For example, on the row i, the entries (i, 5) and (i, 9) can be updated simultaneously,
#' then the entries (i, 1), (i, 3), (i, 8), and so on. 
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{time}{elapsed time to run the MCMC in second.}
#'     \item{iteration}{number of MCMC steps performed.}
#'     \item{posterior}{matrix containing the simulations.}
#'     \item{hyperparms}{return value of `hyperparms`.}
#'     \item{accept.rate}{acceptance rate of zeta.}
#'     \item{G}{last draw of G (row normalized).}
#'     \item{start}{starting values.}
#'     \item{formula}{input value of `formula`.}
#'     \item{contextual}{input value of `contextual`.}
#'     \item{ctrl.mcmc}{return value of `ctrl.mcmc`.}
#' @examples 
#' \dontrun{
#' # Number of groups
#' M             <- 10
#' # size of each group
#' N             <- rep(50,M)
#' # precision parameter for the network formation process
#' lambda        <- 1 
#' 
#' G             <- list()
#' 
#' # individual effects
#' beta          <- c(2,1,1.5) 
#' # contextual effects
#' gamma         <- c(5,-3) 
#' # endogenous effect
#' alpha         <- 0.4
#' # std-dev errors
#' se            <- 1 
#' 
#' prior         <-list()
#' 
#' ## generate network probabilities
#' for (m in 1:M) {
#'   Nm          <- N[m]
#'     c           <- rnorm(Nm*Nm,0,1)
#'     # linking probabilities
#'     Prob        <- matrix(exp(c/lambda)/(1+exp(c/lambda)),Nm) 
#'     # no self-link
#'    diag(Prob)  <- 0 
#'    prior[[m]]  <- Prob
#' }
#' 
#' ## generate data
#' # covariates
#' X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
#' # dependent variable
#' y             <- c()
#' 
#' for (m in 1:M) {
#'   Nm          <- N[m]
#'   # true network
#'   Gm          <- matrix(runif(Nm^2),Nm,Nm) < prior[[m]] 
#'   # no self-link
#'   diag(Gm)    <- rep(0,Nm) 
#'   G[[m]]      <- Gm
#'   rsm         <- rowSums(Gm)
#'   rsm[rsm==0] <- 1
#'   # normalize
#'   Gm          <- Gm/rsm 
#'   # rows index of group m
#'   r2          <- sum(N[1:m])
#'   r1          <- r2 - Nm + 1
#'   # contextual effect
#'   Xm          <- X[r1:r2,]
#'   GXm         <- Gm \%*\% Xm
#'   tmp         <- cbind(rep(1,Nm),Xm)\%*\%beta + GXm\%*\%gamma + rnorm(Nm,0,se)
#'   y[r1:r2]    <- solve(diag(Nm)-alpha*Gm)\%*\%tmp 
#' }
#'   
#' # number of parameters
#' Kv            <- 2*ncol(X) + 1 
#' 
#' # set the hyperparameter
#' # the hyperparameter is a list
#' hyperparms    <- list("dnetwork" = prior) 
#'
#' # launch the MCMC
#' # update entry on A, one by one
#' out1          <- mcmcSAR(y ~ X | X, hyperparms = hyperparms)
#' summary(out1)
#' plot(out1)
#' plot(out1, plot.type = "dens")
#' # update up to 4 entries simultaneously
#' ctrl <- list(block.max = 4, print.level = 2)
#' out2          <- mcmcSAR(y ~ X | X, hyperparms = hyperparms, ctrl.mcmc = ctrl)
#' summary(out2)
#' plot(out2)
#' plot(out2, plot.type = "dens")
#'}
#' @references 
#' Atchadé, Y. F., & Rosenthal, J. S. (2005). On adaptive markov chain monte carlo algorithms. \emph{Bernoulli}, 11(5), 815-828. \url{https://projecteuclid.org/euclid.bj/1130077595}.
#' @references 
#' Boucher, V., & Houndetoungan, A. (2019). Estimating peer effects using partial network data. \emph{Draft avaliable at} \url{https://houndetoungan.wixsite.com/aristide/research}.
#' @importFrom Formula as.Formula
#' @importFrom stats model.frame
#' @seealso \code{\link{sim.IV}}
#' @export
mcmcSAR <- function(formula,
                    contextual,
                    start,
                    hyperparms,
                    G0     = NULL,
                    iteration = 2000L,
                    ctrl.mcmc = list(),
                    data){
  t1           <- Sys.time()
  target       <- ctrl.mcmc$target
  jumpmin      <- ctrl.mcmc$jumpmin
  jumpmax      <- ctrl.mcmc$jumpmax
  print.level  <- ctrl.mcmc$print.level
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
    if (block.max < 1 | block.max > 10) {
      stop("block.max should is less than 1 or greater than 10")
    }
  }
  
  if (is.null(c)) {
    c          <- 0.6
  } else {
    if(length(c) != 1) {
      stop("c in ctrl.mcmc should be a scalar")
    }
  }
  
  if (is.null(print.level)) {
    print.level      <- 1
  } else {
    if(length(print.level) != 1) {
      stop("print.level in ctrl.mcmc should be a scalar")
    }
  }
  
  ctrl.mcmc    <- list(
    target     = target,
    jumpmin    = jumpmin,
    jumpmax    = jumpmax,
    block.max  = block.max,
    print.level= print.level,
    c          = c
  )
  
  if (missing(contextual)) {
    contextual <- FALSE
  }
  
  f.t.data     <- formula.to.data(formula, contextual, data)
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
  ktheta       <- kbeta + kgamma
  
  col.Xone     <- colnames(Xone)
  col.x        <- colnames(X)
  col.post     <- c(col.Xone)
  if (is.null(X)) {
    col.post     <- c(col.post, "Gy", "se2")
  } else {
    col.post     <- c(col.post, paste0("G: ", col.x), "Gy", "se2")
  }
  

  # hyperparameter
  dnetwork     <- hyperparms$dnetwork
  if (is.null(dnetwork)) {
    stop("dnetwork is not defined in hyperparms")
  }
  if (!is.list(dnetwork)) {
    stop("dnetwork in hyperparms is not a list")
  }
  if (length(dnetwork) == 0) {
    stop("dnetwork in hyperparms is an empty list")
  }
  mutheta      <- hyperparms$mutheta
  if (is.null(mutheta)) {
    mutheta    <- rep(0,ktheta)
  }
  invstheta    <- hyperparms$invstheta
  if (is.null(invstheta)) {
    invstheta  <- diag(ktheta)/100
  }
  muzeta       <- hyperparms$muzeta
  if (is.null(muzeta)) {
    muzeta     <- 0
  }
  invszeta     <- hyperparms$invszeta
  if (is.null(invszeta)) {
    invszeta   <- 2
  }
  a            <- hyperparms$a
  if (is.null(a)) {
    a          <- 4.2
  }
  b            <- hyperparms$b
  if (is.null(b)) {
    b          <- (a - 2)/1
  }
  
  hyperparms   <- c(hyperparms, 
                    list(
                      mutheta   = mutheta,
                      invstheta = invstheta,
                      muzeta    = muzeta,
                      invszeta  = invszeta,
                      a         = a,
                      b         = b
                    ))
  
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
    if (missing(start)) {
      start    <- sartpointnoc(G0, M, N, kbeta, y, Xone);
    }
    if (block.max == 1) {
      out      <- peerMCMCnoc(y, Xone, G0, M, N, kbeta, dnetwork, mutheta, invstheta, 
                           muzeta, invszeta, a, b, start, iteration, target, jumpmin, jumpmax,
                           c, print.level)
    } else {
      out      <- peerMCMCblocknoc(y, Xone, G0, M, N, kbeta, dnetwork, mutheta, invstheta, 
                                muzeta, invszeta, a, b, start, iteration, target, jumpmin, jumpmax,
                                c, block.max, print.level)
    }
  } else {
    if (missing(start)) {
      start    <- sartpoint(G0, M, N, kbeta, kgamma, y, X, Xone);
    }
    if (block.max == 1) {
      out      <- peerMCMC(y, X, Xone, G0, M, N, kbeta, kgamma, dnetwork, mutheta, invstheta, 
                            muzeta, invszeta, a, b, start, iteration, target, jumpmin, jumpmax,
                           c, print.level)
    } else {
      out      <- peerMCMCblock(y, X, Xone, G0, M, N, kbeta, kgamma, dnetwork, mutheta, invstheta, 
                                muzeta, invszeta, a, b, start, iteration, target, jumpmin, jumpmax,
                                c, block.max, print.level)
    }
  }
  
  colnames(out$posterior) <- col.post
  
  
  cat("\n\n")
  cat("The program successfully executed \n")
  cat("\n")
  cat("*************SUMMARY************* \n")
  cat("Number of group        : ", M, "\n")
  cat("Iteration              : ", iteration, "\n")
  
  # Print the processing time
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs")) 
  nhours     <- floor(timer/3600)
  nminutes   <- floor((timer-3600*nhours)/60)%%60
  nseconds   <- timer-3600*nhours-60*nminutes
  cat("Elapsed time           : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
  cat("Average acceptance rate: ", out$acceptance, "\n")

  
  out        <- c(list("n.group"     = M,
                       "N"           = c(N),
                       "time"        = timer,
                       "iteration"   = iteration,
                       "posterior"   = out$posterior,
                       "hyperparms"  = hyperparms, 
                       "accept.rate" = out$acceptance, 
                       "G"           = out$G,
                       "start"       = c(start),
                       "formula"     = formula,
                       "contextual"  = contextual,
                       "ctrl.mcmc"   = ctrl.mcmc))
  
  class(out) <- "mcmcSAR"
  return(out)
}