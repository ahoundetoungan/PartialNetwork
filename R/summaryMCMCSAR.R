#' @title Summarizing Bayesian SAR Model 
#' @description Summary and print methods for the class `mcmcSAR`.
#' @param object an object of class "mcmcSAR", output of the function \code{\link{mcmcSAR}}.
#' @param x an object of class "summary.mcmcSAR" or "mcmcSAR, output of the functions \code{\link{summary.mcmcSAR}} and
#' \code{\link{print.summary.mcmcSAR}}.
#' @param alpha (optional, default = 0.95), the significance level of parameter.
#' @param plot.type (optional) a character that indicate if the simulations from the posterior distribution should be printed
#' (if `plot.type = "sim"`) or if the posterior distribution densities should be plotted (`plot.type = "dens"`). The plots can also
#' done using the method \link[graphics]{plot}.
#' @param burnin is the number of MCMC steps which will be considered as burn-in iterations. If `NULL` (default value),
#' the 50% first MCMC steps performed are used as burn-in iterations.
#' @param ... further arguments passed to or from other methods.
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{iteration}{number of MCMC steps performed.}
#'     \item{burnin}{number of MCMC steps which will be considered as burn-in iterations.}
#'     \item{posterior}{matrix (or list of matrices) containing the simulations.}
#'     \item{hyperparms}{return value of `hyperparms`.}
#'     \item{accept.rate}{acceptance rate of zeta.}
#'     \item{prop.net}{proportion of observed network data.}
#'     \item{method.net}{network formation model specification.}
#'     \item{formula}{input value of `formula`.}
#'     \item{alpha}{significance level of parameter.}
#'     \item{ctrl.mcmc}{return value of `ctrl.mcmc`.}
#'     \item{...}{arguments passed to methods.}
#' @details The function is smart and allows all the possible arguments with the functions \link[base]{summary}, 
#' \link[graphics]{plot}, \link[graphics]{par}... such as `col`, `lty`, `mfrow`... \code{summary.mcmcSAR},
#' \code{print.summary.mcmcSAR} and \code{print.mcmcSAR} can be called by \code{summary} or \code{print}.
#' @importFrom stats quantile sd
#' @export

"summary.mcmcSAR" <- function(object, alpha = 0.95, plot.type = NULL, burnin = NULL, ...) {
  stopifnot(inherits(object, "mcmcSAR"))
  pos.theta   <- object$posterior
  pos.rho     <- NULL
  if (object$method.net != "none"){
    if (object$method.net != "latent space") {
      pos.rho <- pos.theta$rho
    }
    pos.theta <- pos.theta$theta
  } 
  
  iteration   <- object$iteration
  stopifnot(alpha > 0 & alpha < 1)
  
  
  if (is.null(burnin)) {
    burnin    <- ceiling(iteration/2)
  }
  
  if (burnin + 2 >= iteration) {
    stop("burnin is too close to the number of MCMC steps performed.")
  }
  
  K              <- ncol(pos.theta)
  pos.theta[,K]  <- sqrt(pos.theta[,K])
  cnames         <- c("Mean", "Std.Error", "Inf CI", "Sup CI", "Sign")
  rnames.t       <- colnames(pos.theta)
  rnames.t[K-1]  <- "Peer effects"
  rnames.t[K]    <- "se"
  spost.t        <- fsum.post(post = pos.theta, burnin = burnin, iteration = iteration, 
                              cnames = cnames, rnames = rnames.t, alpha = alpha) 
  spost.r        <- NULL
  if(!is.null(pos.rho)) {
    rnames.r     <- colnames(pos.rho)
    spost.r      <- fsum.post(post = pos.rho, burnin = burnin, iteration = iteration, 
                              cnames = cnames, rnames = rnames.r, alpha = alpha)
  }
  
  sum.post       <- list("theta" = spost.t, "rho" = spost.r)
  
  out                <- list(n.group    = object$n.group,
                             N           = object$N,
                             iteration   = object$iteration,
                             burnin      = burnin,
                             posterior   = sum.post,
                             hyperparms  = object$hyperparms, 
                             accept.rate = object$accept.rate,
                             prop.net    = object$prop.net,
                             method.net  = object$method.net,
                             formula     = object$formula,
                             alpha       = alpha,
                             ctrl.mcmc   = object$ctrl.mcmc,
                             ...         = ...)
  
  if (!is.null(plot.type)) {
    stopifnot(plot.type %in% c("sim", "dens"))
    object           <- list(x = object, plot.type = plot.type,... = ...)
    tmp.plot         <- list(x = do.call("plot.mcmcSAR", object), ... = ...)
    do.call("print.plot.mcmcSAR", tmp.plot)
  }
  
  class(out)         <- "summary.mcmcSAR"
  out
}

#' @rdname summary.mcmcSAR
#' @export
"print.summary.mcmcSAR"  <- function(x, ...) {
  stopifnot(inherits(x, "summary.mcmcSAR"))
  
  M          <- x$n.group
  N          <- x$N
  sumpost    <- x$posterior$theta
  sumposr    <- x$posterior$rho
  K          <- nrow(sumpost)
  sposttmp   <- c(list(x = sumpost[-K,], ... = ...), x[-(1:12)])
  if (!is.null(sumposr)) {
    sposrtmp <- c(list(x = sumposr, ... = ...), x[-(1:12)])
  }
  sumN       <- sum(N)
  prop.net   <- x$prop.net
  method.net <- x$method.net
  
  
  cat("Bayesian estimation of SAR model", "\n\n")
  cat("Outcome model's formula = ", Reduce(paste, deparse(x$formula$outcome.model)), "\n", sep = "")
  cat("Method: MCMC\n")
  cat("Number of steps performed: ", x$iteration, "\n", sep = "")
  cat("Burn-in: ", x$burnin, "\n", sep = "")
  
  cat("\nPercentage of observed network data: ", prop.net$propG0.obs*100, "%", sep = "")
  cat("\nNetwork formation model: ", method.net, sep = "")
  if (method.net %in% c("logit", "probit")){
    cat("\nFormula = ", Reduce(paste, deparse(x$formula$network.model)), sep = "")
  }
  if (method.net %in% c("latent space")){
    cat("\nPercentage of observed ARD: ", prop.net$propARD*100, "%", sep = "")
  }
  
  cat("\n\nNetwork sampling \nMethod: Gibbs sampler")
  b.max  <- x$ctrl.mcmc$block.max
  r.max  <- ifelse(b.max > 1, paste0("Yes (max size: ",b.max,")") , "No")
  cat("\nUpdate per block: ", r.max, "\n", sep = "")
  if (!is.null(sumposr)) {
    cat("\nNetwork formation model\n")
    do.call("print", sposrtmp)
    rm("sposrtmp")
  }
  cat("\nOutcome model\n")
  do.call("print", sposttmp)
  rm("sposttmp")
  cat("---\nSignificance level: ", round(x$alpha*100),"%\n", sep = "")
  cat("' ' = non signif.  '+' = signif. positive  '-' = signif. negative\n\n")
  
  cat("Error standard-deviation: ", sumpost[K,1], sep = "")
  cat("\nNumber of groups: ", M, sep = "")
  cat("\nTotal sample size: ", sumN, sep = "", "\n\n")
  
  if (method.net == "none") {
    cat("Peer effects acceptance rate: ", x$accept.rate, "\n", sep = "")
  } 
  if (method.net %in% c("probit", "logit")) {
    cat("Peer effects acceptance rate: ", x$accept.rate["alpha"], "\n", sep = "")
    cat("rho acceptance rate         : ", x$accept.rate["rho"], "\n", sep = "")
  }
  if (method.net == "latent space") {
    cat("Peer effects acceptance rate: ", x$accept.rate$alpha, "\n", sep = "")
    cat("rho acceptance rate         : ", mean(x$accept.rate$rho), "\n", sep = "")
  }
  
  invisible(x)
}

#' @rdname summary.mcmcSAR
#' @export
"print.mcmcSAR" <- function(x, ...) {
  stopifnot(inherits(x, "mcmcSAR"))
  print(summary(x, ...))
}

#' @title Plotting estimation of Bayesian SAR model
#' @description Plotting the simulation from the posterior distribution as well as the 
#' density functions of Bayesian SAR model parameter. For more details about the graphical
#' parameter arguments, see \link[graphics]{par}.
#' @param x object of class "mcmcSAR", output of the function \code{\link{mcmcSAR}} or object of class "plot.mcmcSAR", output of the function \code{\link{plot.mcmcSAR}}.
#' @param plot.type character indicating the type of plot: `"sim"` for plotting the simulation
#' from the posterior distribution or `"dens"` for plotting the posterior density functions.
#' @param burnin number of MCMC steps which will be considered as burn-in iterations. If `NULL` (default value),
#' the 50% first MCMC steps performed are used as burn-in iterations.
#' @param which.parms character indicating the parameters whose the posterior distribution will be plotted: `"theta"` for
#' the parameters of the outcome model and `"rho"` for the parameters of the network formation model.
#' @param ... arguments to be passed to methods, such as \link[graphics]{par}.
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{iteration}{number of MCMC steps performed.}
#'     \item{burnin}{number of MCMC steps which will be considered as burn-in iterations.}
#'     \item{posterior}{summary of the posterior distribution to be plotted.}
#'     \item{hyperparms}{return value of `hyperparms`.}
#'     \item{accept.rate}{acceptance rate of zeta.}
#'     \item{propG0.obs}{proportion of observed network data.}
#'     \item{method.net}{network formation model specification.}
#'     \item{formula}{input value of `formula`.}
#'     \item{ctrl.mcmc}{return value of `ctrl.mcmc`.}
#'     \item{which.parms}{return value of `which.parms`.}
#'     \item{plot.type}{type of the plot.}
#'     \item{...}{arguments passed to methods.}
#' @importFrom graphics par
#' @importFrom stats density
#' @export
"plot.mcmcSAR" <- function(x, plot.type = "sim", burnin = NULL, which.parms = "theta", ...) {
  stopifnot(inherits(x, "mcmcSAR"))
  stopifnot(plot.type %in% c("sim", "dens"))
  stopifnot(which.parms %in% c("theta", "rho"))
  if(x$method.net %in% c("none", "latent space") & which.parms == "rho") {
    stop("which.parms = 'rho' is only allowed for logit and probit models")
  }
  
  iteration      <- x$iteration
  
  if (is.null(burnin)) {
    burnin       <- ceiling(iteration/2)
  }
  
  if (burnin + 2 >= iteration) {
    stop("burnin is too close to the number of MCMC steps performed.")
  }
  
  posterior      <- x$posterior
  if (x$method.net != "none"){
    posterior <- posterior[[which.parms]]
  } 
  k              <- ncol(posterior)
  
  x              <- list(n.group     = x$n.group,
                         N           = x$N,
                         iteration   = x$iteration,
                         burnin      = burnin,
                         posterior   = posterior,
                         hyperparms  = x$hyperparms, 
                         accept.rate = x$accept.rate,
                         propG0.obs  = x$propG0.obs,
                         method.net  = x$method.net,
                         formula     = x$formula,
                         ctrl.mcmc   = x$ctrl.mcmc,
                         which.parms = which.parms,
                         plot.type   = plot.type)
  
  if (is.null(match.call()$mfrow)) {
    nrow         <- NULL
    ncol         <- NULL
    if(k == 1){
      nrow       <- 1
      ncol       <- 1
    }
    else{
      if(k == 2){
        nrow     <- 1
        ncol     <- 2
      }
      else{
        if(k > 2 & k < 5){
          nrow   <- 2
          ncol   <- 2
        }
        else{
          nrow   <- ceiling(sqrt(k)); ncol <- ceiling(k/nrow)
        }
      }
    }
    x            <- c(x, list(mfrow = c(nrow, ncol)))
  } 
  
  x              <- c(x, list(... = ...))
  
  class(x)       <- "plot.mcmcSAR"
  x
}


#' @rdname plot.mcmcSAR
#' @export
"print.plot.mcmcSAR" <- function(x, ...) {
  posterior       <- x$posterior
  k               <- ncol(posterior)
  
  burnin          <- x$burnin
  iteration       <- x$iteration
  
  oldpar          <- par(no.readonly = TRUE)   
  on.exit(par(oldpar))    
  
  do.call("par", c(x[-(1:13)], list(...)))
  sim             <- x$plot.type == "sim"
  
  coln            <- colnames(posterior)
  if (x$which.parms == "theta") {
    coln[k - 1]   <- "Peer effects"
    posterior[,k] <- sqrt(posterior[,k])
  } else{
    k             <- k + 1
  }
  
  if (sim) {
    for (i in 1:(k - 1)) {
      tmp         <- c(list(x = posterior[, i],
                            type = "l",
                            main = coln[i],
                            xlab = "",
                            ylab = "",
                            ...  = ...), x[-(1:14)])
      do.call("plot", tmp)
    }
    if (x$which.parms == "theta") {
      tmp         <- c(list(x = posterior[, k],
                            type = "l",
                            main = expression(sigma),
                            xlab = "",
                            ylab = "",
                            ...  = ...), x[-(1:14)])
      do.call("plot", tmp)
    }
  } else {
    posterior     <- posterior[(burnin+1):iteration,,drop = FALSE]
    for (i in 1:(k - 1)) {
      tmp         <- c(list(x    = density(posterior[, i]),
                            main = coln[i],
                            xlab = "",
                            ylab = "",
                            ...  = ...), x[-(1:14)])
      do.call("plot", tmp)
    }
    if (x$which.parms == "theta") {
      tmp         <- c(list(x    = density(posterior[, k]),
                            main = expression(sigma),
                            xlab = "",
                            ylab = "",
                            ...  = ...), x[-(1:13)])
      do.call("plot", tmp)}
  }
  
  
  # par(mfrow = c(1, 1))
  
  invisible(x)
}

fsum      <- function(x, alpha = 0.95) {
  c(mean(x), sd(x), quantile(x, probs = c(0.5*(1 - alpha), 0.5*(1 + alpha))))
} 


percent   <- function(w, digits = 3, format = "f", ...) {
  paste0(formatC(100 * w, format = format, digits = digits, ...), "%")
}

fsum.post <- function(post, burnin, iteration, cnames, rnames, alpha = 0.95) {
  post               <- post[(burnin+1):iteration,, drop = FALSE]
  
  sum.post           <- t(apply(post, 2, function(x) fsum(x, alpha)))
  sum.post           <- as.data.frame(sum.post)
  
  post.sign          <- apply(sum.post, 1, function(x) ifelse(x[3]*x[4] > 0, ifelse(x[1] > 0, "+", "-"), ""))
  sum.post$sign      <- post.sign 
  
  colnames(sum.post) <- cnames
  rownames(sum.post) <- rnames
  sum.post
}
