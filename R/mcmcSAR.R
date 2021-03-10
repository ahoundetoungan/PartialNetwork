#' @title Bayesian Estimator of SAR model
#' @description \code{mcmcSAR} implements the Bayesian estimator of the linear-in-mean SAR model when only the linking probabilities are available or can be estimated.
#' @param formula object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' the `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param  start (optional) vector of starting value of the model parameter as \eqn{(\beta' ~ \gamma' ~ \alpha ~ \sigma^2)'}{(\beta'  \gamma'  \alpha  se2)'},
#' where \eqn{\beta} is the individual variables parameter, \eqn{\gamma} is the contextual variables parameter, \eqn{\alpha} is the peer effect parameter
#' and \eqn{\sigma^2}{se2} the variance of the error term. If the `start` is missing, a Maximum Likelihood estimator will be used, where
#' the network matrix is that given through the argument `G0` (if provided) or generated from it distribution `dnetwork` (see argument `hyperparms`).
#' @param G0.obs list of matrices (or simply matrix if the list contains only one matrix) indicating the part of the network data which is observed. If the (i,j)-th element
#' of the m-th matrix is one, then the element at the same position in the network data will be considered as observed and will not be inferred in the MCMC. In contrast, 
#' if the (i,j)-th element of the m-th matrix is zero, the element at the same position in the network data will be considered as a starting value of the missing true value which will be inferred. 
#' `G0.obs` can also take `none` when no part of the network data is observed (equivalent to the case where all the entries are zeros) and `all` when the network data is fully 
#' observed (equivalent to the case where all the entries are ones).
#' @param G0 list of sub-network matrices (or simply network matrix if there is only one sub-network). `G0` is made up of starting values for the entries with missing network data and observed values for the entries with
#' observed network data. `G0` is optional when `G0.obs = "none"`. 
#' @param mlinks list specifying the network formation model (see details).
#' @param hyperparms (optional) is a list of hyperparameters (see details). 
#' @param ctrl.mcmc list of MCMC controls (see details).
#' @param iteration number of MCMC steps to be performed.
#' @param data optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @details 
#' ## Outcome model
#' The model is given by
#' \deqn{\mathbf{y} = \mathbf{X}\beta + \mathbf{G}\mathbf{X}\gamma + \alpha \mathbf{G}\mathbf{y} + \epsilon.}{y = X\beta + GX\gamma + \alpha Gy + \epsilon.}
#' The parameters to estimate in this model are the matrix \eqn{\mathbf{G}}{G}, the vectors \eqn{\beta}, \eqn{\gamma} and the scalar \eqn{\alpha}, \eqn{\sigma^2}{se2}.
#' Prior distributions are assumed on \eqn{\mathbf{A}}, the adjacency matrix in which \eqn{\mathbf{A}_{ij} = 1}{A[i,j] = 1} if i is  connected to j and
#' \eqn{\mathbf{A}_{ij} = 0}{A[i,j] = 0} otherwise, and on \eqn{\beta}, \eqn{\gamma}, \eqn{\alpha} and \eqn{\sigma^2}{se2}.
#' \deqn{\mathbf{A}_{ij} \sim Bernoulli(\mathbf{P}_{ij})}{A[i,j] ~ Bernoulli(P[i,j])}
#' \deqn{(\beta' ~ \gamma')'|\sigma^2 \sim \mathcal{N}(\mu_{\theta}, \sigma^2\Sigma_{\theta})}{(\beta' ~ \gamma')'|se2 ~ N(mutheta, se2*stheta)}
#' \deqn{\zeta = \log\left(\frac{\alpha}{1 - \alpha}\right) \sim \mathcal{N}(\mu_{\zeta}, \sigma_{\zeta}^2)}{\zeta = log(\alpha/(1 - \alpha)) ~ N(muzeta, szeta)}
#' \deqn{\sigma^2 \sim IG(\frac{a}{2}, \frac{b}{2})}{se2 ~ IG(a/2, b/2)}
#' where \eqn{\mathbf{P}}{P} is the linking probability.
#' ## Network formation model
#' Forthcoming
#' ## Hyperparameters
#' All the hyperparameters can be defined through the argument `hyperparms` (a list) and should be named as follow.
#' \itemize{
#' \item `mutheta`, the prior mean of \eqn{(\beta' ~ \gamma')'|\sigma^2}{(\beta' ~ \gamma')'|se2}. The default value assumes that
#' the prior mean is zero.
#' \item `invstheta` as \eqn{\Sigma_{\theta}^{-1}}{inverse of `stheta`}. The default value is a diagonal matrix with 0.01 on the diagonal.
#' \item `muzeta`, the prior mean of \eqn{\zeta}. The default value is zero.
#' \item `invszeta`, the inverse of the prior variance of \eqn{\zeta} with default value equal to 2.
#' \item `a` and `b` which default values equal to 4.2 and 2.2 respectively. This means for example that the prior mean of \eqn{\sigma^2}{se2} is 1.
#' }
#' Inverses are used for the prior variance through the argument `hyperparms`  in order to allow non informative prior. Set the inverse of the prior
#' variance to 0 is equivalent to assume a non informative prior.
#' ## Adaptive MCMC
#' During the MCMC, the jumping scales of \eqn{\alpha} and \eqn{\rho} are updated following Atchadé and Rosenthal (2005) in order to target the acceptance rate to the `target` value. This
#' requires to set a minimal and a maximal jumping scales through the parameter `ctrl.mcmc`. The parameter `ctrl.mcmc` is a list which can contain the following named components.
#' \itemize{
#' \item{`target`}: the default value is \code{c("alpha" = 0.44, "rho" = 0.27)}.
#' \item{`jumpmin`}: the default value is \code{c("alpha" = 1e-12, "rho" = 1e-12)}. 
#' \item{`jumpmax`}: the default value is \code{c("alpha" = 10, "rho" = 10)}. 
#' \item{`print.level`}: an integer in \{0, 1, 2\} that indicates if the MCMC progression should be printed in the console.
#'  If 0, the MCMC progression is not be printed. If 1 (default value), the progression is printed and if 2,
#'  the simulations from the posterior distribution are printed.
#' \item{`block.max`}: The maximal number of entries that can be updated simultaneously in \eqn{\mathbf{A}}{A}. It might be 
#' more efficient to update simultaneously 2 or 3 entries (see Boucher and Houndetoungan, 2020). 
#' }
#' If `block.max` > 1, several entries are randomly chosen from the same row and updated simultaneously. The number of entries chosen is randomly 
#' chosen between 1 and `block.max`. In addition, the entries are not chosen in order. For example, on the row i, the entries (i, 5) and (i, 9) can be updated simultaneously,
#' then the entries (i, 1), (i, 3), (i, 8), and so on. 
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{time}{elapsed time to run the MCMC in second.}
#'     \item{iteration}{number of MCMC steps performed.}
#'     \item{posterior}{matrix containing the simulations.}
#'     \item{hyperparms}{return value of `hyperparms`.}
#'     \item{mlinks}{return value of `mlinks`.}
#'     \item{accept.rate}{acceptance rate of zeta.}
#'     \item{G}{last draw of G (row normalized).}
#'     \item{start}{starting values.}
#'     \item{formula}{input value of `formula`.}
#'     \item{contextual}{input value of `contextual`.}
#'     \item{ctrl.mcmc}{return value of `ctrl.mcmc`.}
#' @examples 
#' \donttest{
#' #################### EXAMPLE I: WITHOUT NETWORK FORMATION MODEL ####################
#' # Number of groups
#' M             <- 50
#' # size of each group
#' N             <- rep(30,M)
#' # individual effects
#' beta          <- c(2,1,1.5) 
#' # contextual effects
#' gamma         <- c(5,-3) 
#' # endogenous effects
#' alpha         <- 0.4
#' # std-dev errors
#' se            <- 1 
#' # prior distribution
#' prior         <- runif(sum(N*(N-1)))
#' prior         <- vec.to.mat(prior, N, normalise = FALSE)
#' # covariates
#' X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
#' # true network
#' G0            <- sim.network(prior)
#' # normalise 
#' G0norm        <- norm.network(G0)
#' # simulate dependent variable use an external package
#' y             <- CDatanet::simSARnet(~ X, contextual = TRUE, Glist = G0norm, 
#'                                      theta = c(alpha, beta, gamma, se))
#' y             <- y$y
#' # dataset
#' dataset       <- as.data.frame(cbind(y, X1 = X[,1], X2 = X[,2])) 
#' ############### Example I-1: When the network is fully observed ###############
#' out.none1     <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#'                          G0 = G0, data = dataset, iteration = 1e4)
#' summary(out.none1)
#' plot(out.none1)
#' plot(out.none1, plot.type = "dens")
#' 
#' ############### Example I-2: When a part of the network is observed ###############
#' # 60% of the network data is observed
#' G0.obs       <- lapply(N, function(x) matrix(rbinom(x^2, 1, 0.6), x))
#' # replace the non-observed part of the network by 0 (missing links)
#' G0.start     <- lapply(1:M, function(x) G0[[x]]*G0.obs[[x]])
#' # Use network with missing data as the true network
#' out.none2.1  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#'                         G0 = G0.start,   data = dataset, iteration = 1e4)
#' summary(out.none2.1) # the peer effets seem overestimated
#' plot(out.none2.1)
#' plot(out.none2.1, plot.type = "dens")
#' # Infer the missing links in the network data
#' # In this example, the distribution of the network is set fixed
#' out.none2.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = G0.obs,
#'                         G0 = G0.start,   data = dataset, mlinks = list(dnetwork = prior),
#'                         iteration = 2e4)
#' summary(out.none2.2)
#' plot(out.none2.2)
#' plot(out.none2.2, plot.type = "dens")
#' 
#' ############### Example I-3: When only the network distribution is available ###############
#' # Simulate a fictitious network and use as true network
#' G0.tmp       <- sim.network(prior)
#' out.none3.1  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#'                         G0 = G0.tmp, data = dataset, iteration = 1e4)
#' summary(out.none3.1)  # the peer effets seem overestimated
#' plot(out.none3.1)
#' plot(out.none3.1, plot.type = "dens")
#' # Infer the network data
#' out.none3.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "none",
#'                         G0 = G0.tmp, data = dataset, mlinks = list(dnetwork = prior),
#'                         iteration = 2e4)
#' summary(out.none3.2)  
#' plot(out.none3.2)
#' plot(out.none3.2, plot.type = "dens")
#' 
#' 
#' #################### EXAMPLE II: NETWORK FORMATION MODEL: LOGIT ####################
#' # Number of groups
#' M             <- 50
#' # size of each group
#' N             <- rep(30,M)
#' # individual effects
#' beta          <- c(2,1,1.5) 
#' # contextual effects
#' gamma         <- c(5,-3) 
#' # endogenous effects
#' alpha         <- 0.4
#' # std-dev errors
#' se            <- 2 
#' # parameters of the network formation model
#' rho           <- c(-2, -.5, .2)
#' # covariates
#' X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
#' # compute distance between individuals 
#' tmp           <- c(0, cumsum(N))
#' X1l           <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],1])
#' X2l           <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],2])
#' dist.net      <- function(x, y) abs(x - y) # distance function
#' X1.mat        <- lapply(1:M, function(m) {
#'                  matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), N[m])}
#'                  )
#' X2.mat        <- lapply(1:M, function(m) {
#'                  matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), N[m])}
#'                  )
#' # true network
#' covar         <- as.matrix(cbind("Const" = 1, 
#'                                  "dX1"   = mat.to.vec(X1.mat), 
#'                                  "dX2"   = mat.to.vec(X2.mat)))
#' ynet          <- covar %*% rho
#' ynet          <- 1*((ynet + rlogis(length(ynet))) > 0)
#' G0            <- vec.to.mat(ynet, N, normalise = FALSE)
#' G0norm        <- norm.network(G0)
#' # simulate dependent variable use an external package
#' y             <- CDatanet::simSARnet(~ X, contextual = TRUE, Glist = G0norm, 
#'                                      theta = c(alpha, beta, gamma, se))
#' y             <- y$y
#' # dataset
#' dataset       <- as.data.frame(cbind(y, X1 = X[,1], X2 = X[,2])) 
#' ############### Example II-1: When the network is fully observed ###############
#' out.prob1     <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#'                          G0 = G0norm, data = dataset, iteration = 1e4)
#' summary(out.prob1)
#' plot(out.prob1)
#' plot(out.prob1, plot.type = "dens")
#' 
#' ############### Example II-2: When a part of the network is observed ###############
#' # 60% of the network data is observed
#' G0.obs       <- lapply(N, function(x) matrix(rbinom(x^2, 1, 0.6), x))
#' # replace the non-observed part of the network by 0
#' G0.start     <- lapply(1:M, function(x) G0[[x]]*G0.obs[[x]])
#' # Use network with missing data as the true network
#' out.prob2.1  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#'                         G0 = G0.start,   data = dataset, iteration = 1e4)
#' summary(out.prob2.1) # the peer effets seem overestimated
#' plot(out.prob2.1)
#' plot(out.prob2.1, plot.type = "dens")
#' # Infer the missing links in the network data
#' out.prob2.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = G0.obs,
#'                         G0 = G0.start,   data = dataset, 
#'                         mlinks = list(model = "logit", covariates = covar), iteration = 2e4)
#' summary(out.prob2.2)
#' plot(out.prob2.2)
#' plot(out.prob2.2, plot.type = "dens")
#' plot(out.prob2.2, which.parms = "rho")
#' 
#' ############### Example II-3: When only the network distribution is available ###############
#' # Simulate a fictitious network and use as true network
#' G0.tmp       <- sim.network(vec.to.mat(plnorm(covar %*% rho), N))
#' out.prob3.1  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#'                         G0 = G0.tmp, data = dataset, iteration = 1e4)
#' summary(out.prob3.1)  # the peer effets seem overestimated
#' plot(out.prob3.1)
#' plot(out.prob3.1, plot.type = "dens")
#' # Infer the network data
#' # We provide information only about rho
#' Gvec         <- mat.to.vec(G0, ceiled = TRUE)
#' logestim     <- glm(Gvec ~ -1 + covar, family = binomial(link = "logit"))
#' slogestim    <- summary(logestim)
#' estimates    <- list("rho"     = logestim$coefficients, 
#'                      "var.rho" = slogestim$cov.unscaled,
#'                      "N"       = N)
#' out.prob3.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "none", data = dataset, 
#'                         mlinks = list(model = "logit", covariates = covar, estimates = estimates),
#'                         iteration = 2e4)
#' summary(out.prob3.2)  
#' plot(out.prob3.2)
#' plot(out.prob3.2, plot.type = "dens")
#' plot(out.prob3.2, which.parms = "rho")
#' }
#' @importFrom Formula as.Formula
#' @importFrom stats model.frame
#' @importFrom Matrix rankMatrix
#' @importFrom stats glm
#' @importFrom stats pnorm
#' @importFrom stats plogis
#' @importFrom stats binomial
#' @importFrom stats cov
#' @importFrom utils tail
#' @seealso \code{\link{sim.IV}}
#' @export
mcmcSAR <- function(formula,
                    contextual,
                    start,
                    G0.obs,
                    G0         = NULL,
                    mlinks     = list(),
                    hyperparms = list(),
                    ctrl.mcmc  = list(),
                    iteration  = 2000L,
                    data){
  
  t1                <- Sys.time()
  # data
  if (missing(contextual)) {
    contextual <- FALSE
  }
  
  f.t.data     <- formula.to.data(formula, contextual, data)
  formula      <- f.t.data$formula
  Xone         <- f.t.data$Xone
  X            <- f.t.data$X
  y            <- f.t.data$y
  
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
  
  ### model of links
  lmodel       <- toupper(mlinks$model)
  estimates    <- mlinks$estimates
  dZ           <- mlinks$covariates
  mARD         <- mlinks$mARD
  burninARD    <- mlinks$burninARD
  obsARD       <- mlinks$obsARD
  dnetwork     <- mlinks$dnetwork
  murho        <- NULL
  Vrho         <- NULL
  ListIndex    <- NULL
  Krho         <- NULL
  P            <- NULL
  typeprob     <- NULL
  lFdZrho1     <- NULL
  lFdZrho0     <- NULL
  cn.rho       <- NULL
  dest         <- NULL
  zetaest      <- NULL
  if (!is.null(dZ)) {
    dZ         <- as.matrix(dZ)
  }
  tmodel       <- f.tmodel(lmodel, G0.obs, G0, dZ, estimates, dnetwork, obsARD)
  M            <- tmodel$M
  N            <- tmodel$N
  obsARD       <- tmodel$obsARD
  sumG0.obs    <- tmodel$sumG0.obs
  lmodel       <- tmodel$lmodel  #NONE PROBIT LOGIT or LATENT SPACE
  tmodel       <- tmodel$tmodel  #NONE PARTIAL ALL
  
  # target, jumping
  target        <- f.targ.jump(x = ctrl.mcmc$target, sval = c(0.44, 0.27), lmodel) 
  jumpmin       <- f.targ.jump(x = ctrl.mcmc$jumpmin, sval = rep(1e-12, 2), lmodel) 
  jumpmax       <- f.targ.jump(x = ctrl.mcmc$jumpmax, sval = rep(10, 2), lmodel) 
  print.level   <- ctrl.mcmc$print.level
  block.max     <- ctrl.mcmc$block.max
  cpar          <- ctrl.mcmc$c
  if (is.null(cpar)) {
    cpar        <- 0.6
  } 
  if (is.null(block.max)) {
    block.max   <- 1
  }
  if (is.null(print.level)) {
    print.level <- 1
  } 
  block.max    <- block.max[1]
  cpar         <- cpar[1]
  print.level  <- print.level[1]
  block.max    <- round(block.max)
  
  if (block.max < 1 | block.max > 10) {
    stop("block.max should is less than 1 or greater than 10")
  }
  
  ctrl.mcmc    <- list(
    target          = target,
    jumpmin         = jumpmin,
    jumpmax         = jumpmax,
    block.max       = block.max,
    print.level     = print.level,
    c               = cpar
  )
  # hyperparameter
  mutheta      <- hyperparms$mutheta
  invstheta    <- hyperparms$invstheta
  muzeta       <- hyperparms$muzeta
  invszeta     <- hyperparms$invszeta
  a            <- hyperparms$a
  b            <- hyperparms$b
  
  if (is.null(mutheta)) {
    mutheta    <- rep(0,ktheta)
  }
  if (is.null(invstheta)) {
    invstheta  <- diag(ktheta)/100
  }
  if (is.null(muzeta)) {
    muzeta     <- 0
  }
  
  if (is.null(invszeta)) {
    invszeta   <- 2
  }
  if (is.null(a)) {
    a          <- 4.2
  }
  if (is.null(b)) {
    b          <- (a - 2)/1
  }
  
  hyperparms   <- list(
    mutheta   = mutheta,
    invstheta = invstheta,
    muzeta    = muzeta,
    invszeta  = invszeta,
    a         = a,
    b         = b
  )
  
  # model NONE
  if(lmodel == "NONE") {
    if(tmodel == "NONE") {
      if(class(G0.obs) != "list") {
        G0.obs <- lapply(N, function(x) matrix(0, x, x))
      }
    }
    if(tmodel == "ALL") {
      if(class(G0.obs) != "list") {
        G0.obs <- lapply(N, function(x) matrix(1, x, x))
      }
      dnetwork <- G0.obs
    }
    ListIndex  <- fListIndex(dnetwork, G0.obs, M, N)
  }
  
  # model Probit and logit
  if(lmodel %in% c("PROBIT", "LOGIT")) {
    typeprob   <- ifelse(lmodel == "PROBIT", 1, 2)
    Krho       <- ncol(dZ)
    cn.rho     <- colnames(dZ)
    cn.rho     <- unlist(lapply(1:Krho, function(x) ifelse(length(cn.rho[x]) == 0, paste0("covariates.", x), cn.rho[x])))
    
    murho      <- estimates$rho
    Vrho       <- estimates$var.rho
    
    if (is.null(murho) | is.null(Vrho)) {
      G0vec    <- frMceiltoV(G0, N, M)
      Gobsvec  <- as.logical(frMceiltoV(G0.obs, N, M))
      G0vec    <- G0vec[Gobsvec]
      dZtmp    <- dZ[Gobsvec, ]
      myp      <- glm(G0vec ~ -1 + dZtmp, family = binomial(link = ifelse(typeprob == 1, "probit", "logit")))
      smyp     <- summary(myp)
      murho    <- smyp$coefficients[,1]
      Vrho     <- smyp$cov.unscaled
    } 

    pfit       <- ifelse(typeprob == 1, pnorm, plogis)(c(dZ%*%murho))
    
    lFdZrho1   <- log(pfit)
    lFdZrho0   <- log(1- pfit)
    
    dnetwork   <- frVtoM(pfit, N, M)
    
    if(tmodel == "NONE") { #else is necessary partial
      if(class(G0.obs) != "list") {
        G0.obs    <- lapply(N, function(x) matrix(0, x, x))
      }
    }
    ListIndex  <- fListIndex(dnetwork, G0.obs, M, N)
  }  
  
  # latent space model
  if(lmodel == "LATENT SPACE") {
    if (print.level) {
      cat("Estimation of the initial parameters for the latent space model\n")
    }
    murho              <- list()  # mean z nu
    Vrho               <- list()  # covariance z nu
    dnetwork           <- list()
    dest               <- list()
    lzetaest           <- c()
    Krho               <- c()
    P                  <- unique(unlist(lapply(estimates, function(x) x$p)))
    if (length(P) != 1) {
      stop("Latent space models with different hypersphere dimensions")
    }
    
    typeprob            <- unlist(lapply(estimates, function(x){
      if(is.numeric(x$accept.rate$d) & is.numeric(x$accept.rate$zeta)) return(0)
      if(is.character(x$accept.rate$d) & is.numeric(x$accept.rate$zeta)) return(1)
      if(is.numeric(x$accept.rate$d) & is.character(x$accept.rate$zeta)) return(2)
      if(is.character(x$accept.rate$d) & is.character(x$accept.rate$zeta)) return(3)
    }))
    
    for (m in 1:M) {
      T              <- estimates[[m]]$iteration
      z              <- estimates[[m]]$simulations$z
      d              <- estimates[[m]]$simulations$d
      zeta           <- estimates[[m]]$simulations$zeta
      out            <- NULL
      Mst            <- round(T/2) + 1
      if (!is.null(burninARD)){
        Mst          <- burninARD + 1
      }
      if (Mst >= T) {
        stop("BurninARD is too high")
      }
      
      ngraph         <- T - Mst + 1
      if(ngraph < 30) {
        stop("Number of simulations performed in the latent space model (after burnin) is low")
      }
      
      lspaceZNU      <- NULL
      if (is.null(dZ)) {
        # estimate Z and NU
        lspaceZNU      <- flspacerho1(T, P, z, d, zeta, N[m],  Mst)
      } else {
        obsARDm        <- obsARD[[m]]
        iARDm          <- which(obsARDm)
        inonARDm       <- which(!obsARDm)
        XARDm          <- dZ[iARDm,]
        XnonARDm       <- dZ[inonARDm,]
        iARDm          <- iARDm - 1
        inonARDm       <- inonARDm - 1

        N1             <- nrow(XARDm)
        N2             <- nrow(XnonARDm)

        # estimate rho
        lspaceZNU      <- flspacerho2(T, P, z, d, zeta, XARDm, XnonARDm, N1, N2, mARD, Mst, iARDm, inonARDm) 
      }

      # estimates murho and Vrho
      lzeta_z_nu       <- lspaceZNU$rho
      lzeta_z_nu[,1]   <- log(lzeta_z_nu[,1])
      dest[[m]]        <- c(lspaceZNU$degree)
      murhom           <- colMeans(lzeta_z_nu)
      lzetaest[m]      <- murhom[1]
      zm               <- murhom[2:(N[m]*P + 1)]
      num              <- tail(murhom, N[m])
      
      dnetwork[[m]]    <- fdnetARD(exp(lzetaest[m]), zm, num, dest[[m]], N[m], P)
      
      if (typeprob[m] == 0){
        murho[[m]]     <- murhom
        Vrho[[m]]      <- cov(lzeta_z_nu)
        Krho[m]        <- N[m]*(P + 1) + 1
      }
      if (typeprob[m] == 1){
        murho[[m]]     <- murhom[1:(N[m]*P + 1)]
        Vrho[[m]]      <- cov(lzeta_z_nu[,1:(N[m]*P + 1)])
        Krho[m]        <- N[m]*P + 1
      }
      if (typeprob[m] == 2){
        murho[[m]]     <- murhom[-1]
        Vrho[[m]]      <- cov(lzeta_z_nu[,-1])
        Krho[m]        <- N[m]*(P + 1)
      }
      if (typeprob[m] == 3){
        murho[[m]]     <- murhom[2:(N[m]*P + 1)]
        Vrho[[m]]      <- cov(lzeta_z_nu[,2:(N[m]*P + 1)])
        Krho[m]        <- N[m]*P
      }
    }
    
    if(tmodel == "NONE") { #else is necessary partial
      if(class(G0.obs) != "list") {
        G0.obs    <- lapply(N, function(x) matrix(0, x, x))
      }
    }
    ListIndex  <- fListIndex(dnetwork, G0.obs, M, N)
  }
  
  # add value to hyperparms
  if (lmodel != "NONE"){
  hyperparms   <- c(hyperparms, 
                    list(
                     rho     = murho,
                     var.rho = Vrho
                    ))
  }
  
  # Network and starting values
  if (is.null(X)) {
    if (is.null(G0)) {
      tmp      <- flistGnorm1nc(dnetwork, y, Xone, M)
      G0       <- tmp$G
      y        <- tmp$ly
      Xone     <- tmp$lXone
    } else {
      tmp      <- flistGnorm2nc(G0, y, Xone, M)
      G0       <- tmp$G
      y        <- tmp$ly
      Xone     <- tmp$lXone
    }
    if (missing(start)) {
      start    <- sartpointnoc(G0, M, N, kbeta, y, Xone);
    }
  } else {
    if (is.null(G0)) {
      tmp      <- flistGnorm1(dnetwork, y, Xone, X, M)
      G0       <- tmp$G
      y        <- tmp$ly
      Xone     <- tmp$lXone
      X        <- tmp$lX
    } else {
      tmp      <- flistGnorm2(G0, y, Xone, X, M)
      G0       <- tmp$G
      y        <- tmp$ly
      Xone     <- tmp$lXone
      X        <- tmp$lX
    }
    if (missing(start)) {
      start    <- sartpoint(G0, M, N, kbeta, kgamma, y, X, Xone)
    }
  }
  
  start        <- c(start)
  names(start) <- col.post
  
  out          <- NULL
  if (lmodel == "NONE") {
    out        <- SARMCMCnone(y, X, Xone, G0, start, M, N, kbeta, kgamma, dnetwork, ListIndex, mutheta, invstheta, 
                              muzeta, invszeta, a, b, iteration, target, jumpmin, jumpmax,
                              cpar, print.level, block.max)
    colnames(out$posterior)      <- col.post
  }
  if(lmodel %in% c("PROBIT", "LOGIT")) {
    out        <- SARMCMCpl(y, X, Xone, G0, G0.obs, start, M, N, kbeta, kgamma, dnetwork, ListIndex, mutheta, invstheta, 
                            muzeta, invszeta, a, b, dZ, murho, Vrho, Krho, lFdZrho1, lFdZrho0, iteration,
                            target, jumpmin, jumpmax, cpar, print.level, typeprob, block.max)
    colnames(out$posterior$theta) <- col.post
    colnames(out$posterior$rho)   <- cn.rho
  }
  if(lmodel == "LATENT SPACE") {
    out        <- SARMCMCard(y, X, Xone, G0, G0.obs, start, M, N, kbeta, kgamma, dnetwork, ListIndex, mutheta, invstheta, 
                             muzeta, invszeta, a, b, dest, exp(lzetaest), murho, Vrho, Krho, P, iteration,
                             target, jumpmin, jumpmax, cpar, typeprob, print.level, block.max)

    colnames(out$posterior$theta) <- col.post
    out$acceptance$rho            <- c(out$acceptance$rho)

    for (mm in 1:M) {
      outr.nam                          <-  paste0("z", paste0(rep(1:P, each = N[m]), "_", 1:N[m]))
      if (typeprob[mm] %in% c(0, 1)) {
        outr.nam                        <- c("logzeta", outr.nam)
      }
      if (typeprob[mm] %in% c(0, 2)) {
        outr.nam                        <- c(outr.nam,  paste0("nu", 1:N[m]))
      }
      colnames(out$posterior$rho[[mm]]) <- outr.nam
    }

    if(any(out$acceptance$rho == 0)) {
      warning("Network distribution set fixed for at least one sub-network because initial estimate of the distribution is zero for couples who are friends (or one for couples who are not friends)")
    }
  }
  
  
  cat("\n")
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
  
  if (print.level > 0) {
    cat("Elapsed time           : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
    
    if (lmodel == "NONE") {
      cat("Peer effects acceptance rate: ", out$acceptance, "\n", sep = "")
    } 
    if (lmodel %in% c("PROBIT", "LOGIT")) {
      cat("Peer effects acceptance rate: ", out$acceptance["alpha"], "\n", sep = "")
      cat("rho acceptance rate         : ", out$acceptance["rho"], "\n", sep = "")
    }
    if (lmodel == "LATENT SPACE") {
      cat("Peer effects acceptance rate: ", out$acceptance$alpha, "\n", sep = "")
      cat("rho acceptance rate         : ", mean(out$acceptance$rho), "\n", sep = "")
    }
  }
  
  
  out        <- c(list("n.group"     = M,
                       "N"           = c(N),
                       "time"        = timer,
                       "iteration"   = iteration,
                       "posterior"   = out$posterior,
                       "hyperparms"  = hyperparms, 
                       "accept.rate" = out$acceptance, 
                       "propG0.obs"  = sum(sumG0.obs)/sum(N*(N - 1))*100,
                       "method.net"  = tolower(lmodel),
                       "start"       = c(start),
                       "formula"     = formula,
                       "contextual"  = contextual,
                       "ctrl.mcmc"   = ctrl.mcmc))
  
  class(out) <- "mcmcSAR"
  return(out)
}


# MCMC for the model none
SARMCMCnone    <- function(y, X, Xone, G0, start, M, N, kbeta, kgamma, dnetwork, ListIndex, mutheta, invstheta, 
                           muzeta, invszeta, a, b, iteration, target, jumpmin, jumpmax,
                           cpar, print.level, block.max) {
  out          <- NULL
  if (is.null(X)) {
    if (block.max == 1) {
      out      <- peerMCMCnoc_none(y             = y, 
                                   V             = Xone, 
                                   Gnorm         = G0,
                                   prior         = dnetwork, 
                                   ListIndex     = ListIndex, 
                                   M             = M, 
                                   N             = N,
                                   kbeta         = kbeta, 
                                   theta0        = mutheta, 
                                   invsigmatheta = invstheta, 
                                   zeta0         = muzeta,
                                   invsigma2zeta = invszeta, 
                                   a             = a, 
                                   b             = b, 
                                   parms0        = start, 
                                   iteration     = iteration,
                                   target        = target, 
                                   jumpmin       = jumpmin, 
                                   jumpmax       = jumpmax, 
                                   c             = cpar, 
                                   progress      = print.level)
    } else {
      out      <- peerMCMCblocknoc_none(y             = y, 
                                        V             = Xone, 
                                        Gnorm         = G0,
                                        prior         = dnetwork, 
                                        ListIndex     = ListIndex, 
                                        M             = M, 
                                        N             = N,
                                        kbeta         = kbeta, 
                                        theta0        = mutheta, 
                                        invsigmatheta = invstheta, 
                                        zeta0         = muzeta,
                                        invsigma2zeta = invszeta, 
                                        a             = a, 
                                        b             = b, 
                                        parms0        = start, 
                                        iteration     = iteration,
                                        target        = target, 
                                        jumpmin       = jumpmin, 
                                        jumpmax       = jumpmax, 
                                        c             = cpar, 
                                        progress      = print.level,
                                        nupmax        = block.max)
    }
  } else {
    if (block.max == 1) {
      out      <- peerMCMC_none(y             = y,
                                X             = X,
                                Xone          = Xone,
                                Gnorm         = G0,
                                prior         = dnetwork, 
                                ListIndex     = ListIndex, 
                                M             = M, 
                                N             = N,
                                kbeta         = kbeta,
                                kgamma        = kgamma,
                                theta0        = mutheta, 
                                invsigmatheta = invstheta, 
                                zeta0         = muzeta,
                                invsigma2zeta = invszeta, 
                                a             = a, 
                                b             = b, 
                                parms0        = start, 
                                iteration     = iteration,
                                target        = target, 
                                jumpmin       = jumpmin, 
                                jumpmax       = jumpmax, 
                                c             = cpar, 
                                progress      = print.level)
    } else {
      out      <- peerMCMCblock_none(y             = y,
                                     X             = X,
                                     Xone          = Xone,
                                     Gnorm         = G0,
                                     prior         = dnetwork, 
                                     ListIndex     = ListIndex, 
                                     M             = M, 
                                     N             = N,
                                     kbeta         = kbeta,
                                     kgamma        = kgamma,
                                     theta0        = mutheta, 
                                     invsigmatheta = invstheta, 
                                     zeta0         = muzeta,
                                     invsigma2zeta = invszeta, 
                                     a             = a, 
                                     b             = b, 
                                     parms0        = start, 
                                     iteration     = iteration,
                                     target        = target, 
                                     jumpmin       = jumpmin, 
                                     jumpmax       = jumpmax, 
                                     c             = cpar, 
                                     progress      = print.level,
                                     nupmax        = block.max)
    }
  }
  return(out)
}

# MCMC for the models PL
SARMCMCpl      <- function(y, X, Xone, G0, G0.obs, start, M, N, kbeta, kgamma, dnetwork, ListIndex, mutheta, invstheta, 
                           muzeta, invszeta, a, b, dZ, murho, Vrho, Krho, lFdZrho1, lFdZrho0, iteration,
                           target, jumpmin, jumpmax, cpar, print.level,
                           typeprob, block.max) {
  out          <- NULL
  if (is.null(X)) {
    if (block.max == 1) {
      out      <- peerMCMCnoc_pl(y             = y, 
                                 V             = Xone, 
                                 Gnorm         = G0,
                                 G0obs         = G0.obs,
                                 prior         = dnetwork, 
                                 ListIndex     = ListIndex, 
                                 M             = M, 
                                 N             = N,
                                 kbeta         = kbeta, 
                                 theta0        = mutheta, 
                                 invsigmatheta = invstheta, 
                                 zeta0         = muzeta,
                                 invsigma2zeta = invszeta, 
                                 a             = a, 
                                 b             = b, 
                                 dZ            = dZ,
                                 murho         = murho,
                                 Vrho          = Vrho,
                                 Krho          = Krho,
                                 lFdZrho1      = lFdZrho1,
                                 lFdZrho0      = lFdZrho0,
                                 parms0        = start, 
                                 iteration     = iteration,
                                 target        = target,
                                 jumpmin       = jumpmin, 
                                 jumpmax       = jumpmax, 
                                 c             = cpar, 
                                 progress      = print.level,
                                 type          = typeprob)
    } else {
      out      <- peerMCMCblocknoc_pl(y             = y, 
                                      V             = Xone, 
                                      Gnorm         = G0,
                                      G0obs         = G0.obs,
                                      prior         = dnetwork, 
                                      ListIndex     = ListIndex, 
                                      M             = M, 
                                      N             = N,
                                      kbeta         = kbeta, 
                                      theta0        = mutheta, 
                                      invsigmatheta = invstheta, 
                                      zeta0         = muzeta,
                                      invsigma2zeta = invszeta, 
                                      a             = a, 
                                      b             = b, 
                                      dZ            = dZ,
                                      murho         = murho,
                                      Vrho          = Vrho,
                                      Krho          = Krho,
                                      lFdZrho1      = lFdZrho1,
                                      lFdZrho0      = lFdZrho0,
                                      parms0        = start, 
                                      iteration     = iteration,
                                      target        = target,
                                      jumpmin       = jumpmin, 
                                      jumpmax       = jumpmax, 
                                      c             = cpar, 
                                      progress      = print.level,
                                      nupmax        = block.max,
                                      type          = typeprob)
    }
  } else {
    if (block.max == 1) {
      out      <- peerMCMC_pl(y             = y,
                              X             = X,
                              Xone          = Xone,
                              Gnorm         = G0,
                              G0obs         = G0.obs,
                              prior         = dnetwork, 
                              ListIndex     = ListIndex, 
                              M             = M, 
                              N             = N,
                              kbeta         = kbeta,
                              kgamma        = kgamma,
                              theta0        = mutheta, 
                              invsigmatheta = invstheta, 
                              zeta0         = muzeta,
                              invsigma2zeta = invszeta, 
                              a             = a, 
                              b             = b, 
                              dZ            = dZ,
                              murho         = murho,
                              Vrho          = Vrho,
                              Krho          = Krho,
                              lFdZrho1      = lFdZrho1,
                              lFdZrho0      = lFdZrho0,
                              parms0        = start, 
                              iteration     = iteration,
                              target        = target,
                              jumpmin       = jumpmin, 
                              jumpmax       = jumpmax, 
                              c             = cpar,  
                              progress      = print.level,
                              type          = typeprob)
    } else {
      out      <- peerMCMCblock_pl(y             = y,
                                   X             = X,
                                   Xone          = Xone,
                                   Gnorm         = G0,
                                   G0obs         = G0.obs,
                                   prior         = dnetwork, 
                                   ListIndex     = ListIndex, 
                                   M             = M, 
                                   N             = N,
                                   kbeta         = kbeta,
                                   kgamma        = kgamma,
                                   theta0        = mutheta, 
                                   invsigmatheta = invstheta, 
                                   zeta0         = muzeta,
                                   invsigma2zeta = invszeta, 
                                   a             = a, 
                                   b             = b, 
                                   dZ            = dZ,
                                   murho         = murho,
                                   Vrho          = Vrho,
                                   Krho          = Krho,
                                   lFdZrho1      = lFdZrho1,
                                   lFdZrho0      = lFdZrho0,
                                   parms0        = start, 
                                   iteration     = iteration,
                                   target        = target,
                                   jumpmin       = jumpmin, 
                                   jumpmax       = jumpmax, 
                                   c             = cpar, 
                                   progress      = print.level,
                                   nupmax        = block.max,
                                   type          = typeprob)
    }
  }
  return(out)
}

# MCMC for the latent space model
SARMCMCard    <- function(y, X, Xone, G0, G0.obs, start, M, N, kbeta, kgamma, dnetwork, ListIndex, mutheta, invstheta, 
                          muzeta, invszeta, a, b,  dest, zetaest, murho, Vrho, Krho, P, iteration,
                          target, jumpmin, jumpmax, cpar,  typeprob, print.level, block.max) {
  out          <- NULL
  if (is.null(X)) {
    if (block.max == 1) {
      out      <- peerMCMCnoc_ard(y             = y, 
                                  V             = Xone, 
                                  Gnorm         = G0,
                                  G0obs         = G0.obs,
                                  prior         = dnetwork, 
                                  ListIndex     = ListIndex, 
                                  M             = M, 
                                  N             = N,
                                  kbeta         = kbeta, 
                                  theta0        = mutheta, 
                                  invsigmatheta = invstheta, 
                                  zeta0         = muzeta,
                                  invsigma2zeta = invszeta, 
                                  a             = a, 
                                  b             = b, 
                                  d             = dest,
                                  zetaard       = zetaest,   
                                  murho         = murho,
                                  Vrho          = Vrho,
                                  Krho          = Krho,
                                  P             = P,
                                  parms0        = start, 
                                  iteration     = iteration,
                                  target        = target,
                                  jumpmin       = jumpmin, 
                                  jumpmax       = jumpmax, 
                                  c             = cpar, 
                                  type          = typeprob, 
                                  progress      = print.level)
    } else {
      out      <- peerMCMCblocknoc_ard(y             = y, 
                                       V             = Xone, 
                                       Gnorm         = G0,
                                       G0obs         = G0.obs,
                                       prior         = dnetwork, 
                                       ListIndex     = ListIndex, 
                                       M             = M, 
                                       N             = N,
                                       kbeta         = kbeta, 
                                       theta0        = mutheta, 
                                       invsigmatheta = invstheta, 
                                       zeta0         = muzeta,
                                       invsigma2zeta = invszeta, 
                                       a             = a, 
                                       b             = b,
                                       d             = dest,
                                       zetaard       = zetaest,   
                                       murho         = murho,
                                       Vrho          = Vrho,
                                       Krho          = Krho,
                                       P             = P,
                                       parms0        = start, 
                                       iteration     = iteration,
                                       target        = target,
                                       jumpmin       = jumpmin, 
                                       jumpmax       = jumpmax, 
                                       c             = cpar, 
                                       type          = typeprob, 
                                       progress      = print.level,
                                       nupmax        = block.max)
    }
  } else {
    if (block.max == 1) {
      out      <- peerMCMC_ard(y             = y,
                               X             = X,
                               Xone          = Xone,
                               Gnorm         = G0,
                               G0obs         = G0.obs,
                               prior         = dnetwork, 
                               ListIndex     = ListIndex, 
                               M             = M, 
                               N             = N,
                               kbeta         = kbeta,
                               kgamma        = kgamma,
                               theta0        = mutheta, 
                               invsigmatheta = invstheta, 
                               zeta0         = muzeta,
                               invsigma2zeta = invszeta, 
                               a             = a, 
                               b             = b, 
                               d             = dest,
                               zetaard       = zetaest,   
                               murho         = murho,
                               Vrho          = Vrho,
                               Krho          = Krho,
                               P             = P,
                               parms0        = start, 
                               iteration     = iteration,
                               target        = target,
                               jumpmin       = jumpmin, 
                               jumpmax       = jumpmax, 
                               c             = cpar,  
                               type          = typeprob, 
                               progress      = print.level)
    } else {
      out      <- peerMCMCblock_ard(y             = y,
                                    X             = X,
                                    Xone          = Xone,
                                    Gnorm         = G0,
                                    G0obs         = G0.obs,
                                    prior         = dnetwork, 
                                    ListIndex     = ListIndex, 
                                    M             = M, 
                                    N             = N,
                                    kbeta         = kbeta,
                                    kgamma        = kgamma,
                                    theta0        = mutheta, 
                                    invsigmatheta = invstheta, 
                                    zeta0         = muzeta,
                                    invsigma2zeta = invszeta, 
                                    a             = a, 
                                    b             = b, 
                                    d             = dest,
                                    zetaard       = zetaest,   
                                    murho         = murho,
                                    Vrho          = Vrho,
                                    Krho          = Krho,
                                    P             = P,
                                    parms0        = start, 
                                    iteration     = iteration,
                                    target        = target,
                                    jumpmin       = jumpmin, 
                                    jumpmax       = jumpmax, 
                                    c             = cpar, 
                                    type          = typeprob, 
                                    progress      = print.level,
                                    nupmax        = block.max)
    }
  }
  return(out)
}


# This function compute the type of the model and also return N, M
f.tmodel         <- function(lmodel, G0.obs, G0, dZ, estimates, dnetwork, obsARD) {
  # Object's class
  if(!is.null(G0)) {
    if (!is.list(G0)) {
      if (is.matrix(G0)) {
        G0      <- list(G0)
      } else {
        stop("G0 is neither null, nor a matrix nor a list")
      }
    }
  }
  
  if(class(G0.obs) != "list" & class(G0.obs) != "character"){
    if (is.matrix(G0.obs)) {
      G0.obs  <- list(G0.obs)
    } else {
      stop("G0.obs is neither a character, nor a matrix nor a list")
    }
  }
  
  if(class(G0.obs) == "character") {
    G0.obs       <- toupper(G0.obs)
    if(length(G0.obs) != 1 & G0.obs != "ALL" & G0.obs != "NONE"){
      stop("G0.obs should be 'none', 'all', a matrix, or list of matrices")
    }
  }
  
  if(length(lmodel) == 0) {
    lmodel       <- "NONE"
  }
  if(!(lmodel %in% c("PROBIT", "LOGIT", "LATENT SPACE", "NONE"))) {
    stop("model in mlinks should be either none, probit, logit, or latent space")
  }
  
  ## Know the type
  N              <- NULL
  M              <- NULL
  tmodel         <- NULL
  sumG0.obs      <- NULL
  
  # Redefine G0.obs if needed
  if (class(G0.obs) != "character") {
    tmp          <- do.call(rbind, lapply(G0.obs, function(x) c(sum(x > 0) - sum(diag(x) > 0), nrow(x)))) 
    sumG0.obs    <- tmp[,1]   
    N            <- tmp[,2]
    M            <- length(G0.obs)
    tmp          <- all((sumG0.obs - N*(N-1)) == 0) 
    if (all((sumG0.obs == 0))) {
      tmodel     <- "NONE"
    } else {
      if (tmp) {
        tmodel   <- "ALL"
      } else {
        tmodel   <- "PARTIAL"
      }
    }
  } else {
    tmodel       <- G0.obs
  }
  
  
  
  # type of models
  if (tmodel == "ALL") {
    if(is.null(G0)) {
      stop("G0 should be a network matrix or list of sub-network matrices if G0.obs != none")
    }
    M          <- length(G0)
    N          <- unlist(lapply(G0, nrow))
    lmodel     <- "NONE"
    sumG0.obs  <- N*(N-1)
  } else{
    # none
    if (lmodel == "NONE") {
      if(class(dnetwork) != "list") {
        if (is.matrix(dnetwork)) {
          dnetwork  <- list(dnetwork)
        } else {
          stop("dnetwork in mlinks should be a matrix or list of matrices if lmodel = none")
        }
      }
      M        <- length(dnetwork)
      N        <- unlist(lapply(dnetwork, nrow))
    }
    
    # latent space
    if (lmodel == "LATENT SPACE") { 
      if(!is.list(estimates)) {
        stop("For the latent space model, estimates in mlinks should be a list of objects returned by the function mcmcARD")
      }
      if(!all(unlist(lapply(estimates, class)) == "estim.ARD")) {
        stop("for the latent space model, estimates in mlinks should be a list of objects returned by the function mcmcARD")
      }
      
      M        <- length(estimates)
      N        <- unlist(lapply(estimates, function(x)x$n))
      if(!is.null(dZ)){
        if(!is.list(dZ)){
          stop(paste0("for the ",
                      tolower(lmodel),
                      " model, covariates in mlinks should be a list of M covariate matrices (for individuals with and without ARD) to be used to compute distances between individuals."))
        }
        Ntmp     <- unlist(lapply(dZ, nrow))
        if(all(Ntmp <= N)) {
          stop("covariates provided when nrow(covariates[[m]]) is not greater than the ARD sample size")
        }
        if(!is.null(obsARD)) {
          if(!all(unlist(lapply(obsARD, function(x) is.vector(x) & is.logical(x))))){
            stop("one element in obsARD is not a logical vector")
          }
          Ntmp2  <- unlist(lapply(obsARD, length))
          if(any(Ntmp2!=Ntmp)){
            stop("nrow(covariates[[m]]) != length(obsARD[[m]]) for at least one m")
          }
          N      <- Ntmp
        } else {
          N2     <- Ntmp - N
          obsARD <- lapply(1:M, function(x) rep(c(TRUE, FALSE), c(N[x], N2[x])))
        }
      }
    }
    
    # probit and logit
    if (lmodel %in% c("PROBIT", "LOGIT")) { 
      murho    <- estimates$rho
      Vrho     <- estimates$var.rho
      
      if(!is.matrix(dZ)){
        stop(paste0("For the ",
                    tolower(lmodel),
                    " model, covariates in mlinks should be a matrice of covariates of the network formation model."))
      }
      if (tmodel == "NONE") {
        if(is.null(N)) N  <- estimates$N
        if (is.null(murho) | is.null(Vrho) | is.null(N)) {
          stop(paste0("For the ",
                      tolower(lmodel),
                      " model without partial information in G0, estimates in mlinks should be a list containing 'rho', the estimate of rho, 'var.rho', the variance of the estimation, 'M', and 'N', a vector of the number of individuals in each sub-network."))
          
        }
        M                <- length(N)
        sumG0.obs         <- rep(0, M)
      }
      if(nrow(dZ) != sum(N^2 - N)) {
        stop("nrow(covariates) != sum(N^2 - N)")
      }
    }
  }
  
  
  return(list("M" = M, "N" = N, "tmodel" = tmodel, "lmodel" = lmodel, "sumG0.obs" = sumG0.obs, "obsARD" = obsARD))
}

# this function set target and jumping
f.targ.jump      <- function(x, sval, lmodel) {
  if (is.null(x)) {
    if (lmodel == "NONE"){
      x          <- c("alpha" = sval[1])
    } else {
      x          <- c("alpha" = sval[1], "rho" = sval[2])
    }
  } else {
    if (is.null(names(x))) {
      if (lmodel == "NONE"){
        x        <- c("alpha" = x[1])
      } else {
        xtmp     <- x
        x        <- c()
        x[1:2]   <- xtmp   
        names(x) <- c("alpha", "rho")
      }
    } else {
      x           <- x[c("alpha", "rho")]
      x[is.na(x)] <- sval[is.na(x)]
      names(x)    <- c("alpha", "rho")
      if (lmodel == "NONE"){
        x         <- c("alpha" = x[1])
      }
    }
  }
  x
}
