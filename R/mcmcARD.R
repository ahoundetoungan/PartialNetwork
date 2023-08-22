#' @title Estimate network model using ARD
#' @description  \code{mcmcARD} estimates the network model proposed by Breza et al. (2020).
#' @param Y is a matrix of ARD. The entry (i, k) is the number of i's friends having the trait k.
#' @param traitARD is the matrix of traits for individuals with ARD. The entry (i, k) is equal to 1 if i has the trait k and 0 otherwise.
#' @param start is a list containing starting values of `z` (matrix of dimension \eqn{N \times p}), `v` (matrix of dimension \eqn{K \times p}),
#'  `d` (vector of dimension \eqn{N}), `b` (vector of dimension \eqn{K}), `eta` (vector of dimension \eqn{K}) and `zeta` (scalar).
#' @param fixv is a vector setting which location parameters are fixed for identifiability.
#' These fixed positions are used to rotate the latent surface back to a common orientation at each iteration using
#' a Procrustes transformation (see Section Identification in Details).
#' @param consb is a vector of the subset of \eqn{\beta_k}{bk} constrained to the total size (see Section Identification in Details).
#' @param iteration is the number of MCMC steps to be performed.
#' @param sim.d is logical indicating whether the degree `d` will be updated in the MCMC. If `sim.d = FALSE`,
#' the starting value of `d` in the argument `start` is set fixed along the MCMC.
#' @param sim.zeta is logical indicating whether the degree `zeta` will be updated in the MCMC. If `sim.zeta = FALSE`,
#' the starting value of `zeta` in the argument `start` is set fixed along the MCMC.
#' @param hyperparms is an 8-dimensional vector of hyperparameters (in this order) \eqn{\mu_d}{mud},  \eqn{\sigma_d}{sigmad},
#' \eqn{\mu_b}{mub}, \eqn{\sigma_b}{sigmab}, \eqn{\alpha_{\eta}}{alphaeta}, \eqn{\beta_{\eta}}{betaeta}, 
#' \eqn{\alpha_{\zeta}}{alphazeta} and \eqn{\beta_{\zeta}}{betazeta} (see Section Model in Details).
#' @param ctrl.mcmc is a list of MCMC controls (see Section MCMC control in Details).
#' 
#' @details The linking probability is given by
#' ## Model
#' \deqn{P_{ij} \propto \exp(\nu_i + \nu_j + \zeta\mathbf{z}_i\mathbf{z}_j).}{Pij is proportional to exp(nui + nuj + zeta * zi * zj).}
#' McCormick and Zheng (2015) write the likelihood of the model with respect to the spherical coordinate \eqn{\mathbf{z}_i}{zi},
#' the trait locations \eqn{\mathbf{v}_k}{vk}, the degree \eqn{d_i}{di}, the fraction of ties in the network that are
#' made with members of group k \eqn{b_k}{bk}, the trait intensity parameter \eqn{\eta_k}{etak} and \eqn{\zeta}{zeta}. 
#' The following
#' prior distributions are defined.
#' \deqn{\mathbf{z}_i \sim Uniform ~ von ~ Mises-Fisher}{zi ~ Uniform von Mises Fisher}
#' \deqn{\mathbf{v}_k \sim Uniform ~ von ~ Mises-Fisher}{vk ~ Uniform von Mises Fisher}
#' \deqn{d_i \sim log-\mathcal{N}(\mu_d, \sigma_d)}{di ~ log-Normal(mud, sigmad)}
#' \deqn{b_k \sim log-\mathcal{N}(\mu_b, \sigma_b)}{bk ~ log-Normal(mub, sigmab)}
#' \deqn{\eta_k \sim Gamma(\alpha_{\eta}, \beta_{\eta})}{etak ~ Gamma(alphaeta, betaeta)}
#' \deqn{\zeta \sim Gamma(\alpha_{\zeta}, \beta_{\zeta})}{zeta ~ Gamma(alphazeta, betazeta)} 
#' ## Identification
#' For identification, some \eqn{\mathbf{v}_k}{vk} and \eqn{b_k}{bk} need to be exogenously fixed around their given starting value
#' (see McCormick and Zheng, 2015 for more details). The parameter `fixv` can be used
#' to set the desired value for \eqn{\mathbf{v}_k}{vk} while `fixb` can be used to set the desired values for \eqn{b_k}{bk}.\cr
#' ## MCMC control
#' During the MCMC, the jumping scales are updated following Atchade and Rosenthal (2005) in order to target the acceptance rate of each parameter to the `target` values. This
#' requires to set minimal and maximal jumping scales through the parameter `ctrl.mcmc`. The parameter `ctrl.mcmc` is a list which can contain the following named components.
#' \itemize{
#' \item{`target`}: The default value is \code{rep(0.44, 5)}. 
#' The target of every \eqn{\mathbf{z}_i}{zi}, \eqn{d_i}{di}, \eqn{b_k}{bk}, \eqn{\eta_k}{etak} and \eqn{\zeta}{zeta} is  0.44.
#' \item{`jumpmin`}: The default value is \code{c(0,1,1e-7,1e-7,1e-7)*1e-5}. 
#' The minimal jumping of every \eqn{\mathbf{z}_i}{zi} is 0, every \eqn{d_i}{di} is \eqn{10^{-5}}{1e-5}, and every \eqn{b_k}{bk}, \eqn{\eta_k}{etak} and \eqn{\zeta}{zeta} is \eqn{10^{-12}}{1e-12}.
#' \item{`jumpmax`}: The default value is \code{c(100,1,1,1,1)*20}. The maximal jumping scale is 20 except for \eqn{\mathbf{z}_i}{zi} which is set to 2000.
#' \item{`print`}: A logical value which indicates if the MCMC progression should be printed in the console. The default value is `TRUE`.
#' }
#' @return A list consisting of:
#'     \item{n}{dimension of the sample with ARD.}
#'     \item{K}{number of traits.}
#'     \item{p}{hypersphere dimension.}
#'     \item{time}{elapsed time in second.}
#'     \item{iteration}{number of MCMC steps performed.}
#'     \item{simulations}{simulations from the posterior distribution.}
#'     \item{hyperparms}{return value of hyperparameters (updated and non updated).}
#'     \item{accept.rate}{list of acceptance rates.}
#'     \item{start}{starting values.}
#'     \item{ctrl.mcmc}{return value of `ctrl.mcmc`.}
#' @examples 
#' \donttest{
#'   # Sample size
#'   N       <- 500
#'   
#'   # ARD parameters
#'   genzeta <- 1
#'   mu      <- -1.35
#'   sigma   <- 0.37
#'   K       <- 12    # number of traits
#'   P       <- 3     # Sphere dimension
#'   
#'   
#'   # Generate z (spherical coordinates)
#'   genz    <- rvMF(N,rep(0,P))
#'   
#'   # Generate nu  from a Normal distribution with parameters mu and sigma (The gregariousness)
#'   gennu   <- rnorm(N,mu,sigma)
#'   
#'   # compute degrees
#'   gend <- N*exp(gennu)*exp(mu+0.5*sigma^2)*exp(logCpvMF(P,0) - logCpvMF(P,genzeta))
#'   
#'   # Link probabilities
#'   Probabilities <- sim.dnetwork(gennu,gend,genzeta,genz)
#'   
#'   # Adjacency matrix
#'   G <- sim.network(Probabilities)
#'   
#'   # Generate vk, the trait location
#'   genv <- rvMF(K,rep(0,P))
#'   
#'   # set fixed some vk  distant
#'   genv[1,] <- c(1,0,0)
#'   genv[2,] <- c(0,1,0)
#'   genv[3,] <- c(0,0,1)
#'   
#'   # eta, the intensity parameter
#'   geneta   <-abs(rnorm(K,2,1))
#'   
#'   # Build traits matrix
#'   densityatz       <- matrix(0,N,K)
#'   for(k in 1:K){
#'     densityatz[,k] <- dvMF(genz,genv[k,]*geneta[k])
#'   }
#'   
#'   trait       <- matrix(0,N,K)
#'   NK          <- floor(runif(K, 0.8, 0.95)*colSums(densityatz)/apply(densityatz, 2, max))
#'   for (k in 1:K) {
#'     trait[,k]  <- rbinom(N, 1, NK[k]*densityatz[,k]/sum(densityatz[,k]))
#'   }
#'   
#'   # print a percentage of people having a trait
#'   colSums(trait)*100/N
#'   
#'   # Build ARD
#'   ARD         <- G %*% trait
#'   
#'   # generate b
#'   genb        <- numeric(K)
#'   for(k in 1:K){
#'     genb[k]   <- sum(G[,trait[,k]==1])/sum(G)
#'   }
#'   
#'   ############ ARD Posterior distribution ###################
#'   # initialization
#'   d0     <- exp(rnorm(N)); b0 <- exp(rnorm(K)); eta0 <- rep(1,K);
#'   zeta0  <- 05; z0 <- matrix(rvMF(N,rep(0,P)),N); v0 <- matrix(rvMF(K,rep(0,P)),K)
#'   
#'   # We need to fix some of the vk and bk for identification (see Breza et al. (2020) for details).
#'   vfixcolumn      <- 1:6
#'   bfixcolumn      <- c(3, 5)
#'   b0[bfixcolumn]  <- genb[bfixcolumn]
#'   v0[vfixcolumn,] <- genv[vfixcolumn,]
#'   start  <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0, "zeta" = zeta0)
#'   
#'   # MCMC
#'   out   <- mcmcARD(Y = ARD, traitARD = trait, start = start, fixv = vfixcolumn,
#'                    consb = bfixcolumn, iteration = 5000)
#'   
#'   # plot simulations
#'   # plot d
#'   plot(out$simulations$d[,100], type = "l", col = "blue", ylab = "")
#'   abline(h = gend[100], col = "red")
#'   
#'   # plot coordinates of individuals
#'   i <- 123 # individual 123
#'   {
#'     lapply(1:3, function(x) {
#'       plot(out$simulations$z[i, x,] , type = "l", ylab = "", col = "blue", ylim = c(-1, 1))
#'       abline(h = genz[i, x], col = "red")
#'     })
#'   }
#'   
#'   # plot coordinates of traits
#'   k <- 8
#'   {
#'     lapply(1:3, function(x) {
#'       plot(out$simulations$v[k, x,] , type = "l", ylab = "", col = "blue", ylim = c(-1, 1))
#'       abline(h = genv[k, x], col = "red")
#'     })
#'   }}
#' @export
mcmcARD        <- function(Y, traitARD, start, fixv, consb, iteration = 2000L, sim.d = TRUE, sim.zeta = TRUE, hyperparms = NULL, ctrl.mcmc = list()) {
  t1           <- Sys.time()
  
  # hyperparameters
  if (is.null(hyperparms)) {
    hyperparms <- c(0,1,0,1,5,0.5,1,1)
  } else {
    if(length(hyperparms) != 8) {
      stop("hyperparms should be a 8-dimensional vector")
    }
  }
  
  # MCMC control
  target       <- ctrl.mcmc$target
  jumpmin      <- ctrl.mcmc$jumpmin
  jumpmax      <- ctrl.mcmc$jumpmax
  print        <- ctrl.mcmc$print
  c            <- ctrl.mcmc$c
  
  if (is.null(target)) {
    target     <- c(1,1,1,1,1)*0.44
  } else {
    if(length(target) != 5) {
      stop("target in ctrl.mcmc should be a 5-dimensional vector")
    }
  }
  if (is.null(jumpmin)) {
    jumpmin    <- c(0,1,1e-7,1e-7,1e-7)*1e-5
  } else {
    if(length(jumpmin) != 5) {
      stop("jumpmin in ctrl.mcmc should be a 5-dimensional vector")
    }
  }
  if (is.null(jumpmax)) {
    jumpmax    <- c(100,1,1,1,1)*20
  } else {
    if(length(jumpmax) != 5) {
      stop("jumpmax in ctrl.mcmc should be a 5-dimensional vector")
    }
  }
  if (is.null(c)) {
    c          <- 0.6
  } else {
    if(length(c) != 1) {
      stop("c in ctrl.mcmc should be a scalar")
    }
  }
  if (is.null(print)) {
    print     <- TRUE
  } else {
    if(length(print) != 1) {
      stop("print in ctrl.mcmc should be a scalar")
    }
  }
  namjum         <- c("z", "d", "b", "eta", "zeta")
  names(target)  <- namjum
  names(jumpmin) <- namjum
  names(jumpmax) <- namjum
  
  
  ctrl.mcmc    <- list(
    target     = target,
    jumpmin    = jumpmin,
    jumpmax    = jumpmax,
    print      = print,
    c          = c
  )
  
  z0           <- start$z
  v0           <- start$v
  d0           <- start$d
  b0           <- start$b
  eta0         <- start$eta
  zeta0        <- start$zeta
  stopifnot(start$d > 0)
  stopifnot(start$b > 0)
  stopifnot(start$eta > 0)
  stopifnot(start$zeta > 0)
  
  if (is.null(z0) | is.null(v0) | is.null(d0) | is.null(b0) | is.null(eta0) | is.null(zeta0)) {
    stop("Elements in start should be named as z, v, b, d, eta and zeta")
  }
  
  n            <- nrow(Y)
  K            <- ncol(Y)
  p            <- ncol(z0)
  
  out          <- updateGP(Y, traitARD, z0, v0, d0, b0, eta0, zeta0, fixv, consb, iteration, !sim.d, !sim.zeta,
                           hyperparms, target, jumpmin, jumpmax,  c, print)
  
  out$accept.rate$z      <- c(out$accept.rate$z)
  out$accept.rate$d      <- c(out$accept.rate$d)
  out$accept.rate$b      <- c(out$accept.rate$b)
  out$accept.rate$eta    <- c(out$accept.rate$eta)
  
  if(!sim.d) {
    out$accept.rate$d    <- "Fixed"
  }
  if(!sim.zeta) {
    out$accept.rate$zeta <- "Fixed"
  }
  
  zaccept     <- out$accept.rate$z
  daccept     <- out$accept.rate$d
  baccept     <- out$accept.rate$b
  etaaccept   <- out$accept.rate$eta
  zetaaccept  <- out$accept.rate$zeta
  
  colnames(out$hyperparms$updated) <- c("mud", "sigmad", "mub", "sigmab")
  
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs"))
  
  out         <- c(list("n" = n,
                        "K" = K,
                        "p" = p,
                        "time" = timer,
                        "iteration" = iteration),
                   out,
                   list("start" = start, "ctrl.mcmc" = ctrl.mcmc))
  
  class(out)  <- "estim.ARD"
  if(print) {
    cat("\n")
    cat("The program successfully executed \n")
    cat("\n")
    cat("********SUMMARY******** \n")
    cat("n              : ", n, "\n")
    cat("K              : ", K, "\n")
    cat("Dimension      : ", p, "\n")
    cat("Iteration      : ", iteration, "\n")
    
    
    # Print the processing time
    nhours     <- floor(timer/3600)
    nminutes   <- floor((timer-3600*nhours)/60)%%60
    nseconds   <- timer-3600*nhours-60*nminutes
    cat("Elapsed time   : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
    cat("Average acceptance rate \n")
    cat("                      z: ", mean(zaccept), "\n")
    cat("                      d: ", ifelse(sim.d, mean(daccept), daccept), "\n")
    cat("                      b: ", mean(baccept), "\n")
    cat("                    eta: ", mean(etaaccept), "\n")
    cat("                   zeta: ", zetaaccept, "\n")
  }
  
  out
}