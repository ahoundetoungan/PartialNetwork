#' @title Estimate network model using ARD
#' @description  \code{mcmcARD} estimate the network model proposed by McCormick and Zheng (2015).
#' @param Y is a matrix of ARD. The entry (i, k) is the number of i's friends having the trait k.
#' @param traitARD is the matrix of traits for induviduals with ARD. The entry (i, k) is 1 if i has the trait k and 0 otherwise.
#' @param start is a list containing starting values of `z`, `v`, `d`, `b`, `eta` and `zeta`
#' @param fixv is a vector of which location parameters among the \eqn{p} will be set fixed for identifiability.
#' These fixed positions are used to rotate the latent surface back to a common orientation at each iteration using
#' a Procrustes transformation (see McCormick and Zheng, 2015; Breza et al., 2017 and details).
#' @param consb is a vector of the subset of \eqn{\beta_k}{bk} constrained to the total size (see McCormick and Zheng, 2015; Breza et al., 2017 and details).
#' @param iteration is the number of simulation in the MCMC process
#' @param sim.d is logical indicating weather the degree `d` will be updated in the MCMC. If `sim.d = FALSE`,
#' the starting value of `d` in the argument `start` is set fixed along the process.
#' @param sim.zeta s logical indicating weather the degree `zeta` will be updated in the MCMC. If `sim.zeta = FALSE`,
#' the starting value of `zeta` in the argument `start` is set fixed along the process.
#' @param hyperparms is an 8-dimensional vector of hyperparameters containing \eqn{\mu_d}{mud},  \eqn{\sigma_d}{sigmad},
#' \eqn{\mu_b}{mub}, \eqn{\sigma_b}{sigmab}, \eqn{\alpha_{\eta}}{alphaeta}, \eqn{\beta_{\eta}}{betaeta}, 
#' \eqn{\alpha_{\zeta}}{alphazeta} and \eqn{\beta_zeta}{betazeta} (see details).
#' @param target is a 6-dimensional vector of  targeted acceptance rates for  `z`, `v`, `d`, `b`, `eta` and `zeta`. 
#' @param jumpmin is a 6-dimensional vector of the minimal jumping scales for  `z`, `v`, `d`, `b`, `eta` and `zeta`. 
#' @param jumpmax is a 6-dimensional vector of the miximal jumping scales for  `z`, `v`, `d`, `b`, `eta` and `zeta`. 
#' @param progress is logical; if TRUE, the progression will be printed in the console.
#' @details ---- VINCENT could you plz describe breafly the story behind the model (before we put the formula)? And any usefull information should be in the details.
#' For example why should we fixed set some traits location and some bk and how to set them. -----
#' \deqn{P_{ij} \propto \nu_i + \nu_j + \zeta\mathbf{z}_i\mathbf{z}_j}{Pij propto nui + nuj + zeta * zi * zj}
#' McCormick and Zheng (2015) write the likelihood of the model with respect to the spherical coordinate \eqn{\mathbf{z}_i}{zi},
#' the trait locations \eqn{\mathbf{v}_k}{vk}, the degree \eqn{d_i}{di}, the fraction of ties in the network that are
#' made with members of group k \eqn{b_k}{bk}, the trait intensity parameter \eqn{\eta_k}{etak} and \eqn{\zeta}{zeta}. The following
#' prior distributions are defined.
#' \deqn{\mathbf{z}_i \sim \text{Uniform ~ von Mises FIsher}}{zi ~ Uniform von Mises Fisher}
#' \deqn{\mathbf{v}_k \sim \text{Uniform ~ von Mises FIsher}}{vk ~ Uniform von Mises Fisher}
#' \deqn{d_i \sim log-\mathcal{N}(\mu_d, \sigma_d)}{di ~ log-Normal(mud, sigmad)}
#' \deqn{b_k \sim log-\mathcal{N}(\mu_b, \sigma_b)}{bk ~ log-Normal(mub, sigmab)}
#' \deqn{\eta_k \sim Gamma(\alpha_{\eta}, \beta_{\eta})}{etak ~ Gamma(alphaeta, betaeta)}
#' \deqn{\zeta \sim Gamma(\alpha_{\zeta}, \beta_{\zeta})}{zeta ~ Gamma(alphazeta, betazeta)}
#' During the MCMC, the jumping scales are updated following Atchadé and Rosenthal (2005) in order to target the acceptance rate of each parameter to `target`.
#' @return A list consisting of:
#'     \item{n}{dimension of the sample with ARD}
#'     \item{K}{number of traits}
#'     \item{p}{hypersphere dimension}
#'     \item{iteration}{number of iterations in the MCMC}
#'     \item{time}{Elapsed time in second}
#'     \item{simulations}{simulations from the posterior distribution}
#'     \item{hyperparameters}{vector of hyperparameters}
#'     \item{Acceptance.rate}{list of acceptance rate}
#' @examples 
#' \donotrun{
#' EXAMPLE HERE, BUT AS THE EXAMPLE IS LONG WHITH HOW WE GENERATE THE DATA, 
#' WE COULD ALSO DIRECT USERS TO THE READEME PAGE ONLINE WHICH GIVE MORE DETAILS. WHAT DO YOU THINK (VINCENT)} 
#' @references Atchadé, Y. F., & Rosenthal, J. S. (2005). On adaptive markov chain monte carlo algorithms. Bernoulli, 11(5), 815-828.
#' @references Breza, E., Chandrasekhar, A. G., McCormick, T. H., & Pan, M. (2017). Using aggregated relational data to feasibly
#'  identify network structure without network data (No. w23491). National Bureau of Economic Research.
#' @references McCormick, T. H., & Zheng, T. (2015). Latent surface models for networks using Aggregated Relational Data. 
#' Journal of the American Statistical Association, 110(512), 1684-1695.
#' @export
mcmcARD        <- function(Y, traitARD, start, fixv, consb, iteration = 2000, sim.d = TRUE, sim.zeta = TRUE, hyperparms = NULL, target = NULL, jumpmin = NULL, jumpmax = NULL, progress = TRUE) {
  t1           <- Sys.time()
  if (is.null(hyperparms)) {
    hyperparms <- c(0,1,0,1,5,0.5,1,1)
  }
  if (is.null(target)) {
    target     <- c(1,1,1,1,1)*0.44
  }
  if (is.null(jumpmin)) {
    jumpmin    <- c(1,1,1,1,1)*1e-12
  }
  if (is.null(jumpmax)) {
    jumpmax    <- c(10,1,1,1,1)*20
  }
  
  c            <- 0.6
  
  
  z0           <- start$z
  v0           <- start$v
  d0           <- start$d
  b0           <- start$b
  eta0         <- start$eta
  zeta0        <- start$zeta
  
  if (is.null(z0) | is.null(v0) | is.null(d0) | is.null(b0) | is.null(eta0) | is.null(zeta0)) {
    stop("Elements in start should be named as z, v, b, d, eta and zeta")
  }
    
  n            <- nrow(Y)
  p            <- ncol(z0)
    
  out          <- updateGP(Y, traitARD, z0, v0, d0, b0, eta0, zeta0, fixv, consb, iteration, !sim.d, !sim.zeta,
                           hyperparms, target, jumpmin, jumpmax,  c, progress)
  
  zaccept     <- out$`Acceptance rate`$z
  daccept     <- out$`Acceptance rate`$d
  baccept     <- out$`Acceptance rate`$b
  etaaccept   <- out$`Acceptance rate`$eta
  zetaaccept  <- out$`Acceptance rate`$zeta
  
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs"))
  
  out         <- c(list("n" = n,
                         "K" = K,
                         "p" = p,
                         "iteration" = iteration,
                         "time" = timer),
                    out)
  class(out)  <- "estim.ARD"
  cat("\n\n")
  cat("The program successfully executed \n")
  cat("\n")
  cat("********SUMMARY******** \n")
  cat("n              : ", n, "\n")
  cat("K              : ", k, "\n")
  cat("Dimension      : ", p, "\n")
  cat("Iteration      : ", iteration, "\n")
  
  # Print the processing time
  nhours     <- floor(timer/3600)
  nminutes   <- floor((timer-3600*nhours)/60)%%60
  nseconds   <- timer-3600*nhours-60*nminutes
  cat("Elapsed time   : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
  cat("Average acceptance rate \n")
  cat("                      z: ", mean(zaccept), "\n")
  cat("                      d: ", mean(daccept), "\n")
  cat("                      b: ", mean(baccept), "\n")
  cat("                    eta: ", mean(etaaccept), "\n")
  cat("                   zeta: ", zetaaccept, "\n")
  
  out
}