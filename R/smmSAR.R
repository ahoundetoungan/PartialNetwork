#' @title Simulated Method of Moments (SMM) Estimator of SAR model
#' @description \code{smmSAR} implements the Simulated Method of Moments (SMM) estimator of the linear-in-mean SAR model when only the linking probabilities are available or can be estimated.
#' @param formula object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | gy | gx1 + gx2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables, `gy` is the average of `y` among friends, and
#' `gx1`, `gx2` are the contextual observed variables. If `gy` is observed and `gx1`, `gx2` are not, the formula should be
#' \code{y ~ x1 + x2 | gy}. If `gy` is not observed and `gx1`, `gx2` are, the formula should be \code{y ~ x1 + x2 || gx1 + gx2}. If `gy`, `gx1`, and `gx2` are not observed, the 
#' the formula should simply be \code{y ~ x1 + x2}.
#' @param  contextual logical; if true, this means that all individual variables will be set as contextual variables. In contrast \code{\link{mcmcSAR}},
#' `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is not equivalent to set formula as `y ~ x1 + x2 || gx1 + gx2`. `formula = y ~ x1 + x2` means that `gy`, `gx1`, and `gx2` 
#' are not observed and `contextual = TRUE` means that the estimated model includes contextual effects.
#' @param fixed.effects logical; if true, group heterogeneity is included as fixed effects.
#' @param dnetwork a list, where the m-th elements is the matrix of link probability in the m-th sub-network. 
#' @param W is the weighted-matrix in the objective function of the SMM.
#' @param smm.ctr is the list of some control parameters (see details). 
#' @param cond.var logical; if true the estimator variance conditional on `dnetwork` will be computed.
#' @param data optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If missing, the variables are taken from \code{environment(formula)}, typically the environment from which `smmSAR` is called.
#' @details 
#' The parameter `smm.ctr` is the list of some control parameters such as:
#'  * `R` numbers of draws R (in the package, we assume S = 1 and T = 1);
#'  * `iv.power` number of powers of the network matrix `G` to be used to construct instruments;
#'  * `opt.tol` optimization tolerance that will be used in \link[stats]{optimize};
#'  * `smoother` (logical) which indicates if draws should be performed using the smoother simulator;
#'  * `h` bandwith of the smoother (required if `smoother = TRUE`);
#'  * `print` (logical) indicates if the optimization process should be printed step by step.
#'
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{time}{elapsed time to run the SMM in second.}
#'     \item{estimates}{vector of estimated parameters.}
#'     \item{formula}{input value of `formula`.}
#'     \item{contextual}{input value of `contextual`.}
#'     \item{fixed.effects}{input value of `fixed.effects`.}
#'     \item{smm.ctr}{input value of `smm.ctr`.}
#'     \item{details}{other details of the model.}
#' @examples 
#' # Number of groups
#' M        <- 10
#' # size of each group
#' N        <- rep(20,M)
#' # covariates
#' X        <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
#' # network formation model parameter
#' rho      <- c(-0.8, 0.2, -0.1)
#' # individual effects
#' beta     <- c(2, 1, 1.5, 5, -3)
#' # endogenous effects
#' alpha    <- 0.4
#' # std-dev errors
#' se       <- 1
#' # network
#' tmp      <- c(0, cumsum(N))
#' X1l      <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],1])
#' X2l      <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],2])
#' dist.net <- function(x, y) abs(x - y)
#' X1.mat   <- lapply(1:M, function(m) {
#'   matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), N[m])})
#' X2.mat   <- lapply(1:M, function(m) {
#'   matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), N[m])})
#' Xnet     <- as.matrix(cbind("Const" = 1,
#'                             "dX1"   = mat.to.vec(X1.mat),
#'                             "dX2"   = mat.to.vec(X2.mat)))
#' ynet     <- Xnet %*% rho
#' ynet     <- c(1*((ynet + rlogis(length(ynet))) > 0))
#' G0       <- vec.to.mat(ynet, N, normalise = FALSE)
#' # normalise
#' G0norm   <- norm.network(G0)
#' # Matrix GX
#' GX       <- peer.avg(G0norm, X)
#' # simulate the dependent variable use an external package
#' y        <- CDatanet::simsar(~ X + GX, Glist = G0norm,
#'                              theta = c(alpha, beta, se))
#' Gy       <- y$Gy
#' y        <- y$y
#' # build dataset
#' dataset           <- as.data.frame(cbind(y, X, Gy, GX))
#' colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2")
#' nNet      <- nrow(Xnet) # network formation model sample size
#' Aobs      <- sample(1:nNet, round(0.3*nNet)) # We observed 30%
#' # We can estimate rho using the gml function from the stats package
#' logestim  <- glm(ynet[Aobs] ~ -1 + Xnet[Aobs,], family = binomial(link = "logit"))
#' slogestim <- summary(logestim)
#' rho.est   <- logestim$coefficients
#' rho.var   <- slogestim$cov.unscaled # we also need the covariance of the estimator
#' 
#' d.logit     <- lapply(1:M, function(x) {
#'   out       <- 1/(1 + exp(-rho.est[1] - rho.est[2]*X1.mat[[x]] -
#'                             rho.est[3]*X2.mat[[x]]))
#'   diag(out) <- 0
#'   out})
#' smm.logit   <- smmSAR(y ~ X1 + X2, dnetwork = d.logit, contextual = TRUE,
#'                       smm.ctr  = list(R = 100L, print = TRUE), data = dataset)
#' summary(smm.logit, dnetwork = d.logit, data = dataset)
#' @export
smmSAR <- function(formula,
                   contextual    = FALSE,
                   fixed.effects = FALSE,
                   dnetwork,
                   W             = "identity",
                   smm.ctr       = list(R = 30L, iv.power = 2L, opt.tol  = 1e-4, smoother = FALSE, print = FALSE),
                   cond.var      = TRUE,
                   data){
  
  t1           <- Sys.time()
  sysSeed      <- .GlobalEnv$.Random.seed
  if(is.null(sysSeed)) set.seed(0)
  seed         <- .GlobalEnv$.Random.seed 
  on.exit({
    if (is.null(sysSeed)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      .GlobalEnv$.Random.seed <- sysSeed 
    }
  })
  
  if(!inherits(W, "character")) W <- as.matrix(W)
  stopifnot(inherits(W, c("character", "matrix", "array")))
  
  # controls
  nR           <- smm.ctr$R; if(is.null(nR)) nR <- 30L
  # nS           <- smm.ctr$S
  # nT           <- smm.ctr$T
  iv.power     <- smm.ctr$iv.power; if(is.null(iv.power)) iv.power <- 2L
  opt.tol      <- smm.ctr$opt.tol; if(is.null(opt.tol)) opt.tol <- 1e-4
  smoother     <- smm.ctr$smoother; if(is.null(smoother)) smoother <- FALSE
  hN           <- smm.ctr$h
  wprint       <- smm.ctr$print; if(is.null(wprint)) wprint <- FALSE
  nS           <- 1
  nT           <- 1
  
  stopifnot(iv.power >= 1)
  if(contextual){
    if(iv.power < 2) stop("iv.power should be at least 2 if the model includes contextual effects")
  }
  stopifnot(inherits(smoother, "logical"))
  if(smoother){
    if(is.null(hN)) stop("h is required for the smoother simulator")
  } else{
    hN         <- 0 #unused value to avoid having hN = NULL (error in c++)
  }
  # data
  env.f        <- environment(formula)
  f.t.data     <- formula.to.data.smm(formula = formula, data = data, fixed.effects = fixed.effects) 
  formula      <- f.t.data$formula; environment(formula) <- env.f
  X1           <- f.t.data$X
  y            <- f.t.data$y
  GX2          <- f.t.data$GX
  Gy           <- f.t.data$Gy
  GXobs        <- !is.null(GX2)
  Gyobs        <- !is.null(Gy)
  col.x1       <- colnames(X1)
  col.gx2      <- colnames(GX2)
  intercept    <- ("(Intercept)" %in% col.x1)
  Kx1          <- length(col.x1) 
  if(!GXobs){
    col.gx2    <- paste0("G: ", col.x1[(1 + intercept):Kx1])
  }
  Kx2          <- length(col.gx2)  
  if((Kx1 - intercept) != Kx2) stop("The number of observed contextual variables does not suit")
  X2           <- X1[,(1 + intercept):Kx1, drop = FALSE]
  
  #sizes
  N            <- sapply(dnetwork, nrow)
  M            <- length(N)
  Nsum         <- sum(N)
  Ncum         <- c(0, cumsum(N))
  Ilist        <- lapply(N, diag)
  Pm           <- iv.power - 1
  ninstr       <- Kx1 + iv.power*Kx2
  
  # weight matrix
  if(inherits(W, "character")){
    W          <- diag(ninstr)
  }
  
  # Other
  ctr.b     <- list(R = nR, S = nS, T = nT, distr = dnetwork, Ilist = Ilist, y = y, X1 = X1, X2 = X2, 
                    W = W, smoother = smoother, hN = hN, Kx1 = Kx1, Kx2 = Kx2, M = M, N = N, Pm = Pm, Ncum = Ncum)
  ctr       <- c(ctr.b, list(seed = seed, lower = -1, upper = 1, tol = opt.tol, ninstr = ninstr))
  fbeta     <- NULL
  c.sol     <- FALSE
  fmvzH     <- NULL
  
  # Type = 0
  if(GXobs & Gyobs){
    nS         <- NULL
    nT         <- NULL
    c.sol      <- TRUE  #close solution
    if(fixed.effects){
      if(contextual){
        V      <- as.matrix(cbind(X1, GX2))
        ctr.b  <- c(list(Gy = Gy, GX2 = GX2, V = V, ninstr = ninstr), ctr.b)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("S", "T", "Ilist", "X1", "X2"))]
        fbeta  <- falbeta0fe
        fmvzH  <- fmvzetaH0fe
      } else{
        ctr.b  <- c(list(Gy = Gy, GX2 = GX2, V = X1, ninstr = ninstr), ctr.b)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("S", "T", "Ilist", "X1", "X2", "Kx2"))]
        fbeta  <- falbeta0ncfe
        fmvzH  <- fmvzetaH0ncfe
      }
    } else{
      if(contextual){
        ctr.b  <- c(list(Gy = Gy, GX2 = GX2, V = as.matrix(cbind(X1, GX2)), ninstr = ninstr), ctr.b)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("S", "T", "Ilist", "X1", "X2"))]
        fbeta  <- falbeta0
        fmvzH  <- fmvzetaH0
      } else{
        ctr.b  <- c(list(Gy = Gy, GX2 = GX2, V = X1, ninstr = ninstr), ctr.b)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("S", "T", "Ilist", "X1", "X2", "Kx2"))]
        fbeta  <- falbeta0nc
        fmvzH  <- fmvzetaH0nc
      }
    }
  }
  
  # Type = 1
  if(GXobs & !Gyobs){
    if(fixed.effects){
      if(contextual){
        V      <- as.matrix(cbind(X1, GX2))
        ctr.b  <- c(list(V = V), ctr.b)
        fopt   <- ifelse(wprint, optim1fepr, optim1fe)
        ctr    <- c(list(f = fopt, V = V), ctr)
        fbeta  <- fbeta1fe
        fmvzH  <- fmvzetaH1fe
      } else {
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("Kx2"))]
        fopt   <- ifelse(wprint, optim1ncfepr, optim1ncfe)
        ctr    <- ctr[!(names(ctr) %in% c("Kx2"))]
        ctr    <- c(list(f = fopt), ctr)
        fbeta  <- fbeta1ncfe
        fmvzH  <- fmvzetaH1ncfe
      }
    } else {
      if(contextual){
        V      <- as.matrix(cbind(X1, GX2))
        ctr.b  <- c(list(V = V), ctr.b)
        fopt   <- ifelse(wprint, optim1pr, optim1)
        ctr    <- c(list(f = fopt, V = V), ctr)
        fbeta  <- fbeta1
        fmvzH  <- fmvzetaH1
      } else {
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("Kx2"))]
        fopt   <- ifelse(wprint, optim1ncpr, optim1nc)
        ctr    <- ctr[!(names(ctr) %in% c("Kx2"))]
        ctr    <- c(list(f = fopt), ctr)
        fbeta  <- fbeta1nc
        fmvzH  <- fmvzetaH1nc
      }
    }
  }
  
  # Type = 2
  if(!GXobs & Gyobs){
    c.sol      <- TRUE  #close solution
    if(fixed.effects){
      if(contextual){
        nT     <- NULL
        ctr.b  <- c(list(Gy = Gy, ninstr = ninstr), ctr.b)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("T","Ilist"))]
        fbeta  <- falbeta2fe
        fmvzH  <- fmvzetaH2fe
      } else {
        nS     <- NULL
        nT     <- NULL
        ctr.b  <- c(list(Gy = Gy, ninstr = ninstr), ctr.b)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("S", "T", "Kx2", "Ilist"))]
        fbeta  <- falbeta2ncfe
        fmvzH  <- fmvzetaH2ncfe
      }
    } else {
      if(contextual){
        nT     <- NULL
        ctr.b  <- c(list(Gy = Gy, ninstr = ninstr), ctr.b)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("T","Ilist"))]
        fbeta  <- falbeta2
        fmvzH  <- fmvzetaH2
      } else {
        nS     <- NULL
        nT     <- NULL
        ctr.b  <- c(list(Gy = Gy, ninstr = ninstr), ctr.b)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("S", "T", "Kx2", "Ilist"))]
        fbeta  <- falbeta2nc
        fmvzH  <- fmvzetaH2nc
      }
    }
  }
  
  # Type = 3
  if(!GXobs & !Gyobs){
    if(fixed.effects){
      if(contextual){
        fopt   <- ifelse(wprint, optim3fepr, optim3fe)
        ctr    <- c(list(f = fopt), ctr)
        fbeta  <- fbeta3fe
        fmvzH  <- fmvzetaH3fe
      } else{
        fopt   <-ifelse(wprint, optim3ncfepr, optim3ncfe)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("Kx2"))]
        ctr    <- ctr[!(names(ctr) %in% c("Kx2"))]
        ctr    <- c(list(f = fopt), ctr)
        fbeta  <- fbeta3ncfe
        fmvzH  <- fmvzetaH3ncfe
      }
    } else{
      if(contextual){
        fopt   <- ifelse(wprint, optim3pr, optim3)
        ctr    <- c(list(f = fopt), ctr)
        fbeta  <- fbeta3
        fmvzH  <- fmvzetaH3
      } else{
        fopt   <- ifelse(wprint, optim3ncpr, optim3nc)
        ctr.b  <- ctr.b[!(names(ctr.b) %in% c("Kx2"))]
        ctr    <- ctr[!(names(ctr) %in% c("Kx2"))]
        ctr    <- c(list(f = fopt), ctr)
        fbeta  <- fbeta3nc
        fmvzH  <- fmvzetaH3nc
      }
    }
  }
  
  opt          <- NULL
  theta        <- NULL
  if(c.sol){
    theta      <- do.call(get("fbeta"), ctr.b)
  } else{
    opt        <- do.call(get("optimize"), ctr)
    ealpha     <- opt$minimum
    assign(".Random.seed", seed, envir = .GlobalEnv) #Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function
    Day        <- rep(0, ninstr)
    Ra         <- matrix(0, ninstr, Kx1 + Kx2*contextual)
    tmp        <- c(list(alpha = ealpha, Day = Day, Ra = Ra), ctr.b)
    theta      <- c(ealpha, do.call(get("fbeta"), tmp))
  }
  
  tname        <- c("Gy", col.x1)
  theta        <- c(theta)
  if(contextual) tname <- c(tname, col.gx2)
  names(theta) <- tname
  
  # conditional variance
  avgrad       <- NULL
  avm          <- NULL
  avmm         <- NULL
  if(cond.var){
    ctr.b      <- c(list(alpha = theta[1], beta = theta[-1]), ctr.b)
    if(!("ninstr" %in% names(ctr.b))){
      ctr.b    <- c(list(ninstr = ninstr), ctr.b)
    }
    assign(".Random.seed", seed, envir = .GlobalEnv)#Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function
    tmp        <- do.call(fmvzH, ctr.b)
    avgrad     <- tmp$derM
    avm        <- tmp$sumM/Nsum
    avmm       <- tmp$sumMM/Nsum
  }
  
  # ctr
  if(!smoother) hN           <- NULL
  smm.ctr                    <- list(R = nR)
  # if(!is.null(nS)) smm.ctr$S <- nS
  # if(!is.null(nT)) smm.ctr$T <- nT
  if(!c.sol) smm.ctr$opt.tol <- opt.tol
  smm.ctr$print              <- wprint
  smm.ctr$iv.power           <- iv.power
  smm.ctr$smoother           <- smoother
  if(smoother) smm.ctr$h     <- hN
  
  # Print the processing time
  t2           <- Sys.time()
  timer        <- as.numeric(difftime(t2, t1, units = "secs")) 
  nhours       <- floor(timer/3600)
  nminutes     <- floor((timer-3600*nhours)/60)%%60
  nseconds     <- timer-3600*nhours-60*nminutes
  
  # details
  infos        <- list(iv.power     = iv.power,
                       W            = W,
                       smoother     = smoother,
                       h            = hN,
                       optimize     = opt,
                       GXobs        = GXobs, 
                       Gyobs        = Gyobs, 
                       av.grad.m    = avgrad, 
                       av.m         = c(avm), 
                       "av.m%*%t(m)"= avmm,
                       time         = c("Elapsed times (seconds)" = timer),
                       seed         = seed)
  
  if (wprint > 0) {
    cat("Elapsed time: ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
  }
  out          <- list(n.group       = M,
                       N             = N,
                       estimates     = theta,
                       formula       = formula,
                       contextual    = contextual,
                       fixed.effects = fixed.effects,
                       smm.ctr       = smm.ctr,
                       details       = infos)
  class(out)   <-"smmSAR"
  out
}

#' @title Summarizing SMM Estimation of SAR model
#' @description Summary and print methods for the class `smmSAR`.
#' @param object an object of class "smmSAR", output of the function \code{\link{smmSAR}}.
#' @param x an object of class "summary.smmSAR" or "smmSAR", output of the functions \code{\link{summary.smmSAR}} or
#' \code{\link{smmSAR}}.
#' @param .fun,.args are used to simulate from the distribution of `dnetwork`. `.fun` is the simulator function
#' where `.args` is a list of its arguments. Typically `do.call(.fun, .args)` is supposed to simulate one `dnetwork` from
#' the distribution.
#' @param sim the number of simulations of `dnetwork`.
#' @param ncores the number of cores to be used for the simulation. Use a lot of cores for fast simulations.
#' @param dnetwork a list, where the m-th elements is the matrix of link probability in the m-th sub-network. 
#' @param data optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If missing, the variables are taken from \code{environment(formula)}, typically the environment from which `smmSAR` is called.
#' @param ... further arguments passed to or from other methods.
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{estimates}{vector of estimated parameters.}
#'     \item{formula}{input value of `formula`.}
#'     \item{contextual}{input value of `contextual`.}
#'     \item{fixed.effects}{input value of `fixed.effects`.}
#'     \item{smm.ctr}{input value of `smm.ctr`.}
#'     \item{details}{other details of the model.}
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doRNG "%dorng%"
#' @importFrom stats var
#' @export
"summary.smmSAR" <- function(object, 
                             .fun, 
                             .args, 
                             sim    = 30,
                             ncores = 1,
                             dnetwork, 
                             data,
                             ...){
  stopifnot(inherits(object, "smmSAR"))
  details       <- object$details
  derM          <- details$av.grad.m
  aveMM         <- details$`av.m%*%t(m)`
  aveM          <- details$av.m
  cond.var      <- !(is.null(derM)|is.null(aveMM)|is.null(aveM))
  W             <- details$W
  if(!cond.var | !missing(.fun)) {
    if(missing(dnetwork)|missing(data)) {
      stop("`dnetwork` and `data` should be provided to the `summary` method if variance components were not computed in `smmSAR`.")
    }
  }
  if(!missing(.args) & missing(.fun)){
    stop("`.args` is defined while `.fun` is missing.")
  }
  
  N             <- object$N
  M             <- object$n.group
  if(M < 2) stop("Inference is not possible with one group")
  Nsum          <- sum(N)

  seed          <- details$seed
  sysSeed       <- .GlobalEnv$.Random.seed
  on.exit({
    if (is.null(sysSeed)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      .GlobalEnv$.Random.seed <- sysSeed 
    }
  })
  
  lmethod       <- !cond.var | !missing(.fun)
  varcov        <- NULL
  if(!cond.var | !missing(.fun)) {
    iv.power      <- details$iv.power
    smoother      <- details$smoother
    hN            <- details$h 
    if(smoother){
      if(is.null(hN)) stop("h is required for the smoother simulator")
    } else{
      hN          <- 0 #unused value to avoid having hN = NULL (error in c++)
    }
    formula       <- object$formula
    fixed.effects <- object$fixed.effects
    contextual    <- object$contextual
    
    # data
    f.t.data     <- formula.to.data.smm(formula = formula, data = data, fixed.effects = fixed.effects) 
    formula      <- f.t.data$formula
    X1           <- f.t.data$X
    y            <- f.t.data$y
    GX2          <- f.t.data$GX
    Gy           <- f.t.data$Gy
    GXobs        <- !is.null(GX2)
    Gyobs        <- !is.null(Gy)
    col.x1       <- colnames(X1)
    col.gx2      <- colnames(GX2)
    intercept    <- ("(Intercept)" %in% col.x1)
    Kx1          <- length(col.x1) 
    if(!GXobs){
      col.gx2    <- paste0("G: ", col.x1[(1 + intercept):Kx1])
    }
    Kx2          <- length(col.gx2)  
    if((Kx1 - intercept) != Kx2) stop("The number of observed contextual variables does not suit")
    X2           <- X1[,(1 + intercept):Kx1, drop = FALSE]
    
    # controls
    nR           <- object$smm.ctr$R
    nS           <- 1#object$smm.ctr$S
    nT           <- 1#object$smm.ctr$T
    
    #sizes
    N            <- sapply(dnetwork, nrow)
    M            <- length(N)
    if(M < 2) stop("Inference is not possible with one group")
    Nsum         <- sum(N)
    Ncum         <- c(0, cumsum(N))
    Ilist        <- lapply(N, diag)
    Pm           <- iv.power - 1
    ninstr       <- Kx1 + iv.power*Kx2
    
    # Parameters
    alpha        <- object$estimates[1]
    beta         <- object$estimates[-1]
    
    # Functions fmvzeta and fmvzetaH
    fmvzeta      <- NULL
    fmvzetaH     <- NULL
    Afmvzeta     <- list(alpha = alpha, beta = beta, R = nR, distr = dnetwork, y = y, W = W, smoother = smoother, 
                         hN = hN, Kx1 = Kx1, ninstr = ninstr, M = M, N = N, Pm = Pm, Ncum = Ncum)
    # Type = 0
    if(GXobs & Gyobs){
      if(fixed.effects){
        if(contextual){
          fmvzeta   <- fmvzeta0fe
          fmvzetaH  <- fmvzetaH0fe
          Afmvzeta  <- c(Afmvzeta, list(Gy = Gy, GX2 = GX2, V = as.matrix(cbind(X1, GX2)), Kx2 = Kx2))
        } else{
          fmvzeta   <- fmvzeta0ncfe
          fmvzetaH  <- fmvzetaH0ncfe
          Afmvzeta  <- c(Afmvzeta, list(Gy = Gy, GX2 = GX2, V = X1))
        }
      } else{
        if(contextual){
          fmvzeta   <- fmvzeta0
          fmvzetaH  <- fmvzetaH0
          Afmvzeta  <- c(Afmvzeta, list(Gy = Gy, GX2 = GX2, V = as.matrix(cbind(X1, GX2)), Kx2 = Kx2))
        } else{
          fmvzeta   <- fmvzeta0nc
          fmvzetaH  <- fmvzetaH0nc
          Afmvzeta  <- c(Afmvzeta, list(Gy = Gy, GX2 = GX2, V = X1))
        }
      }
    }
    # Type = 1
    if(GXobs & !Gyobs){
      if(fixed.effects){
        if(contextual){
          fmvzeta   <- fmvzeta1fe
          fmvzetaH  <- fmvzetaH1fe
          Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2, 
                                        V = as.matrix(cbind(X1, GX2)), Kx2 = Kx2))
        } else {
          fmvzeta   <- fmvzeta1ncfe
          fmvzetaH  <- fmvzetaH1ncfe
          Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2))
        }
      } else {
        if(contextual){
          fmvzeta   <- fmvzeta1
          fmvzetaH  <- fmvzetaH1
          Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2, 
                                        V = as.matrix(cbind(X1, GX2)), Kx2 = Kx2))
        } else {
          fmvzeta   <- fmvzeta1nc
          fmvzetaH  <- fmvzetaH1nc
          Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2))
        }
      }
    }
    # Type = 2
    if(!GXobs & Gyobs){
      if(fixed.effects){
        if(contextual){
          fmvzeta   <- fmvzeta2fe
          fmvzetaH  <- fmvzetaH2fe
          Afmvzeta  <- c(Afmvzeta, list(S = nS, X1 = X1, X2 = X2, Gy = Gy, Kx2 = Kx2))
        } else {
          fmvzeta   <- fmvzeta2ncfe
          fmvzetaH  <- fmvzetaH2ncfe
          Afmvzeta  <- c(Afmvzeta, list(X1 = X1, X2 = X2, Gy = Gy))
        }
      } else {
        if(contextual){
          fmvzeta   <- fmvzeta2
          fmvzetaH  <- fmvzetaH2
          Afmvzeta  <- c(Afmvzeta, list(S = nS, X1 = X1, X2 = X2, Gy = Gy, Kx2 = Kx2))
        } else {
          fmvzeta   <- fmvzeta2nc
          fmvzetaH  <- fmvzetaH2nc
          Afmvzeta  <- c(Afmvzeta, list(X1 = X1, X2 = X2, Gy = Gy))
        }
      }
    }
    # Type = 3
    if(!GXobs & !Gyobs){
      if(fixed.effects){
        if(contextual){
          fmvzeta   <- fmvzeta3fe
          fmvzetaH  <- fmvzetaH3fe
          Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2, Kx2 = Kx2))
        } else{
          fmvzeta   <- fmvzeta3ncfe
          fmvzetaH  <- fmvzetaH3ncfe
          Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2))
        }
      } else{
        if(contextual){
          fmvzeta   <- fmvzeta3
          fmvzetaH  <- fmvzetaH3
          Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2, Kx2 = Kx2))
        } else{
          fmvzeta   <- fmvzeta3nc
          fmvzetaH  <- fmvzetaH3nc
          Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2))
        }
      }
    }
    
    # Estimate H0 and OMEGA
    H0        <- NULL
    if(cond.var){
      tderM   <- t(derM)
      H0      <- solve(tderM %*% W %*% derM, tderM)
    } else{
      on.exit({
        if (is.null(sysSeed)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else {
          .GlobalEnv$.Random.seed <- sysSeed 
        }
      })
      assign(".Random.seed", seed, envir = .GlobalEnv)#Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function
      tmp     <- do.call(fmvzetaH, Afmvzeta)
      derM    <- tmp$derM
      tderM   <- t(derM)
      aveMM   <- tmp$sumMM/Nsum
      aveM    <- tmp$sumM/Nsum
      H0      <- solve(tderM %*% W %*% derM, tderM)
      details <- c(details, list(av.grad.m = derM, av.m = c(aveM), "av.m%*%t(m)"= aveMM))
    }
    
    OMEGA   <- NULL
    if(missing(.fun)) {
      OMEGA <- aveMM - Nsum*aveM %*% t(aveM)/M
    } else{
      on.exit({
        if (is.null(sysSeed)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else {
          .GlobalEnv$.Random.seed <- sysSeed 
        }
      })
      stopifnot(sim > 5)
      # Construct cluster
      cl    <- makeCluster(ncores)
      
      # After the function is run, close the cluster.
      on.exit(stopCluster(cl))
      
      # Register parallel backend
      registerDoParallel(cl)
      # assign(".Random.seed", seed, envir = .GlobalEnv)#Note that I restore the system seed (if it exists) using on.exit, see the beginning of the function
      tmp   <- foreach(i = 1:sim, .packages  = "PartialNetwork") %dorng% {
        fOMEGA(.fun, .args, fmvzeta, Afmvzeta, M)}
      
      VZ    <- Reduce('+', lapply(1:sim, function(x) tmp[[x]]$VZ))/(sim*Nsum)
      EZ    <- t(sapply(1:sim, function(x) tmp[[x]]$EZ))
      OMEGA <- VZ + var(EZ)/Nsum
    } 
    varcov  <- H0 %*% W %*% OMEGA %*% t(W) %*% t(H0)/Nsum
  } else {
    tderM   <- t(derM)
    H0      <- solve(tderM %*% W %*% derM, tderM)
    OMEGA   <- aveMM - Nsum*aveM %*% t(aveM)/M
    varcov  <- H0 %*% W %*% OMEGA %*% t(W) %*% t(H0)/Nsum
  }
  
  
  out       <- list(n.group       = object$n.group,
                    N             = object$N,
                    estimates     = object$estimates,
                    cov           = varcov,
                    formula       = object$formula,
                    contextual    = object$contextual,
                    fixed.effects = object$fixed.effects,
                    smm.ctr       = object$smm.ctr,
                    details       = details,
                    ...           = ...)
  class(out) <- "summary.smmSAR"
  out
}


fOMEGA <- function(.fun, .args, fmvzeta, Afmvzeta, M) {
  Afmvzeta$distr  <- do.call(.fun, .args)
  tmp   <- do.call(fmvzeta, Afmvzeta)
  sumMM <- tmp$sumMM
  sumM  <- tmp$sumM
  list(VZ = sumMM - sumM %*% t(sumM)/M, EZ = sumM)
}

#' @rdname summary.smmSAR
#' @export
"print.summary.smmSAR"  <- function(x, ...) {
  stopifnot(inherits(x, "summary.smmSAR"))
  
  M          <- x$n.group
  N          <- x$N
  estimates  <- x$estimates
  sumN       <- sum(N)
  tmp        <- fcoefficients(x$estimates, sqrt(diag(x$cov)))
  out_print  <- tmp$out_print
  out        <- tmp$out
  out_print  <- c(list(out_print), x[-(1:9)], list(...))
  
  nR         <- x$smm.ctr$R
  nS         <- 1#x$smm.ctr$S
  nT         <- 1#x$smm.ctr$T
  
  cat("Simulated Method of Moments estimation of SAR model", "\n\n")
  cat("Formula = ", Reduce(paste, deparse(x$formula)), "\n\n", sep = "")
  cat("Contextual effects: ", ifelse(x$contextual, "Yes", "No"), "\n", sep = "")
  cat("Fixed effects: ", ifelse(x$fixed.effects, "Yes", "No"), "\n\n", sep = "")
  cat("Network details\n")
  cat("GX ", ifelse(x$details$GXobs, "Observed", "Not Observed"), "\n", sep = "")
  cat("Gy ", ifelse(x$details$Gyobs, "Observed", "Not Observed"), "\n", sep = "")
  cat("Number of groups: ", M, "\n", sep = "")
  cat("Sample size     : ", sum(N), "\n\n", sep = "")

  
  #cat("Simulation tuning parameters:\n")
  #cat("R = ", nR, ifelse(!is.null(nS), paste0(", S = ", nS), ""), ifelse(!is.null(nT), paste0(", T = ", nT), "")  , "\n\n", sep = "")
  cat("Simulation settings\n")
  cat("R = ", nR, "\n", sep = "")
  cat("Smoother : ", x$details$smoother, sep = "")
  if(x$details$smoother)cat(" (h = ", x$details$h, ")", sep = "")
  cat("\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  out        <- list(n.group       = x$n.group,
                     N             = x$N,
                     estimates     = x$estimates,
                     cov           = x$cov,
                     coefficients  = out,
                     formula       = x$formula,
                     contextual    = x$contextual,
                     fixed.effects = x$fixed.effects,
                     smm.ctr       = x$smm.ctr)
  class(out) <- "print.summary.smmSAR"
  invisible(out)
}

#' @rdname summary.smmSAR
#' @export
"print.smmSAR"  <- function(x,
                            dnetwork, 
                            .fun, 
                            .args, 
                            sim    = NULL,
                            ncores = 1,
                            data, 
                            ...) {
  stopifnot(inherits(x, "smmSAR"))
  print(summary(object   = x,
                dnetwork = dnetwork, 
                .fun     = .fun, 
                .args    = .args, 
                sim      = sim,
                ncores   = ncores,
                data     = data, ...))
}