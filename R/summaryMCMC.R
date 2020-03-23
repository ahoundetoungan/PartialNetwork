#' @title Summarizing Bayesian SAR Model 
#' @description summary method and print method for class `mcmcSAR`.
#' @param object an object of class "mcmcSAR", output of the function \code{\link{mcmcSAR}}.
#' @param x an object of class "summary.mcmcSAR" or "mcmcSAR, output of thes functions \code{\link{summary.mcmcSAR}} and
#' \code{\link{print.summary.mcmcSAR}}.
#' @param alpha (optional, defaul = 0.95), the significance level of parameter.
#' @param plot.type (optional) a character that indicate if the simulations from the posterior distribution should be printed
#' (if `plot.type = "sim"`) or if the posterior distribution densities should be plotted (`plot.type = "dens"`). The plots can also
#' done using the method \link[graphics]{plot}.
#' @param burnin is the number of MCMC steps which will be considered as burn-in iterations. If `NULL` (default value),
#' the 50\% first MCMC steps performed are used as burn-in iterations.
#' @param ... further arguments passed to or from other methods.
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{iteration}{number of MCMC steps performed.}
#'     \item{burnin}{number of MCMC steps which will be considered as burn-in iterations.}
#'     \item{posterior}{matrix containing the simulations.}
#'     \item{hyperparms}{return value of `hyperparms`.}
#'     \item{accept.rate}{acceptance rate of zeta.}
#'     \item{formula}{input value of `formula`.}
#'     \item{alpha}{significance level of parameter.}
#'     \item{ctrl.mcmc}{return value of `ctrl.mcmc`.}
#'     \item{...}{arguments passed to methods.}
#' @examples 
#' \dontrun{
#' # Let us consider the example in mcmcSAR function
#' out1          <- mcmcSAR(y ~ X | X, hyperparms = hyperparms, ctrl.mcmc = ctrl)
#' # Print the summary
#' summary(out)
#' # Print summary with plot and change significance level 
#' summary(out, alpha = 0.90, plot.type = "sim", col = "blue")
#' # Print summary with plot and change significance level and the layout of the plot
#' summary(out, alpha = 0.90, plot.type = "sim", col = "blue", mfrow = c(4, 4))
#' # All the possible argument with the function summary and the function plot can be passed
#' }
#' @details The function is smart and allows all the possible arguments with the functions \link[base]{summary}, 
#' \link[graphics]{plot}, \link[graphics]{par}... such as `col`, `lty`, `mfrow`... \code{summary.mcmcSAR},
#' \code{print.summary.mcmcSAR} and \code{print.mcmcSAR} can be called by \code{summary} or \code{print}.
#' @importFrom stats quantile sd
#' @export

"summary.mcmcSAR" <- function(object, alpha = 0.95, plot.type = NULL, burnin = NULL, ...) {
  if (!(class(object) == "mcmcSAR")) {
    stop("object is not an mcmcSAR class object")
  }
  posterior <- object$posterior
  iteration <- object$iteration
  stopifnot(alpha > 0 & alpha < 1)
  
  
  if (is.null(burnin)) {
    burnin    <- ceiling(iteration/2)
  }
  
  if (burnin + 2 >= iteration) {
    stop("burnin is too close to the number of MCMC steps performed.")
  }
  posterior <- posterior[(burnin+1):iteration,, drop = FALSE]
  
  K                  <- ncol(posterior)
  posterior[,K]      <- sqrt(posterior[,K])
  
  sum.post  <- t(apply(posterior, 2, function(x) fsum(x, alpha)))
  sum.post  <- as.data.frame(sum.post)
  
  post.sign      <- apply(sum.post, 1, function(x) ifelse(x[3]*x[4] > 0, ifelse(x[1] > 0, "+", "-"), ""))
  sum.post$sign <- post.sign 
  
  colnames(sum.post) <- c("Mean", "Std.Error", "Inf CI", "Sup CI", "Sign")
  rname              <- rownames(sum.post)
  rname[K - 1]       <- "Peer effect"
  rname[K]           <- "se"
  rownames(sum.post) <- rname
  
  out                <- list(n.group     = object$n.group,
                            N           = object$N,
                            iteration   = object$iteration,
                            burnin      = burnin,
                            posterior   = sum.post,
                            hyperparms  = object$hyperparms, 
                            accept.rate = object$accept.rate,
                            formula     = object$formula,
                            alpha       = alpha,
                            ctrl.mcmc   = object$ctrl.mcmc,
                            ...         = ...)
  
  if (!is.null(plot.type)) {
    stopifnot(plot.type %in% c("sim", "dens"))
    object           <- list(object = object, ... = ...)
    tmp.plot         <- list(x = do.call("plot.mcmcSAR", object), ... = ...)
    do.call("print.plot.mcmcSAR", tmp.plot)
  }
  
  class(out)         <- "summary.mcmcSAR"
  out
}

#' @rdname summary.mcmcSAR
#' @export
"print.summary.mcmcSAR"  <- function(x, ...) {
  stopifnot(class(x) == "summary.mcmcSAR")
  
  M          <- x$n.group
  N          <- x$N
  sumpost    <- x$posterior
  K          <- nrow(sumpost)
  sumposttmp <- c(list(x = sumpost[-K,], ... = ...), x[-(1:10)])
  rname      <- rownames(sumpost)
  intercept  <- "(Intercept)" %in% rname
  sumN   <- sum(N)
  dnsum <- colSums(do.call("rbind", lapply(1:M, function(m){
    d <- x$hyperparms$dnetwork[[m]]
    c(sum(d==0) - N[m],
      sum(d==1),
      sum(d<0.20) - N[m],
      sum(d<0.10) - N[m],
      sum(d<0.05) - N[m],
      sum(d<0.01) - N[m],
      sum(d>0.80),
      sum(d>0.90),
      sum(d>0.95),
      sum(d>0.99))
  })))/sum(N*(N - 1))
  
  dnsum1 <- c(rep("<", 4), "=")
  dnsum2 <- c(0.20, 0.10, 0.05, 0.01, 0.00)
  dnsum3 <- dnsum[c(3, 4, 5, 6, 1)]
  dnsum4 <- c(rep(">", 4), "=")
  dnsum5 <- c(0.80, 0.90, 0.95, 0.99, 1.00)
  dnsum6 <- dnsum[c(7, 8, 9, 10, 2)]
  dnsum0 <- rep("", 5)

  
  dnsum  <- data.frame(dnsum1, dnsum2, percent(dnsum3), 
                 dnsum0, dnsum0, dnsum0, dnsum0, dnsum0, dnsum0,
                 dnsum0, dnsum0, dnsum0, dnsum0, dnsum0, dnsum0,
                 dnsum4, dnsum5, percent(dnsum6))
  names(dnsum) <- c(rep("", 2), "proportion", rep("", 14), "proportion")
  
  dnsum  <- c(list(x = dnsum, row.names = FALSE, ... = ...), x[-(1:10)])
  cat("Bayesian estimation of SAR model", "\n\n")
  cat("Formula = ", Reduce(paste, deparse(x$formula)), "\n")
  cat("Method: MCMC\n")
  cat("Number of steps performed: ", x$iteration, "\n")
  cat("Burn-in: ", x$burnin, "\n\n")
  
  cat("Linking probabilities\n")
  do.call("print", dnsum)
  rm(dnsum)
  cat("\nNetwork Sampling \nMethod: Gibbs sampler")
  b.max  <- x$ctrl.mcmc$block.max
  r.max  <- ifelse(b.max > 1, paste0("Yes (max size: ",b.max,")") , "NO")
  cat("\nUpdate per block:", r.max, "\n")
  
  cat("\nPosterior distribution summary \n")
  do.call("print", sumposttmp)
  rm("sumposttmp")
  cat("---\nSignificance level: ", round(x$alpha*100),"%\n")
  cat("( )non signif.  (+)signif. positive  (-)signif. negative\n\n")
  
  cat("Error standard-deviation: ", sumpost[K,1])
  cat("\nNumber of groups:", M)
  cat("\nTotal sample size:", sumN)
  cat("\nPeer effect acceptance rate:", x$accept.rate, "\n")
  
  
    
  invisible(x)
}

#' @rdname summary.mcmcSAR
#' @export
"print.mcmcSAR" <- function(x, ...) {
  stopifnot(class(x) == "mcmcSAR")
  print(summary(x, ...))
}

#' @title Plotting esitmation of Bayesian SAR model
#' @description Plotting the simulation from the posterior distribution as well as the 
#' density functions of Baysian SAR model parameter. For more details about the graphical
#' parameter arguments, see \link[graphics]{par}.
#' @param x an object of class "mcmcSAR", output of the function \code{\link{mcmcSAR}}.
#' @param x an object of class "plot.mcmcSAR", output of thes functions \code{\link{plot.mcmcSAR}}.
#' @param plot.type a character that indicate the type of plot: "sim" for plotting the simulation
#' from the posterior distribution or "dens" for plotting the posterior density functions.
#' @param burnin is the number of MCMC steps which will be considered as burn-in iterations. If `NULL` (default value),
#' the 50\% first MCMC steps performed are used as burn-in iterations.
#' @param ... arguments to be passed to methods, such as \link[graphics]{par}.
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{iteration}{number of MCMC steps performed.}
#'     \item{burnin}{number of MCMC steps which will be considered as burn-in iterations.}
#'     \item{posterior}{summary of the posterior distribution.}
#'     \item{hyperparms}{return value of `hyperparms`.}
#'     \item{accept.rate}{acceptance rate of zeta.}
#'     \item{formula}{input value of `formula`.}
#'     \item{plot.type}{type of the plot.}
#'     \item{ctrl.mcmc}{return value of `ctrl.mcmc`.}
#'     \item{...}{arguments passed to methods.}
#' @examples 
#' \dontrun{
#' # Let us consider the example in mcmcSAR function
#' out1          <- mcmcSAR(y ~ X | X, hyperparms = hyperparms, ctrl.mcmc = ctrl)
#' # Print the summary
#' plot(out)
#' # Print summary with plot and change significance level 
#' plot(out, plot.type = "den", col = "blue")
#' # Print summary with plot and change significance level and the layout of the plot
#' plot(out, plot.type = "dens", col = "blue", mfrow = c(4, 4))
#' # All the possible argument with the function summary and the function plot can be passed
#' }
#' @importFrom graphics par
#' @export
"plot.mcmcSAR" <- function(x, plot.type = "sim", burnin = NULL, ...) {
  stopifnot(class(x) == "mcmcSAR")
  stopifnot(plot.type %in% c("sim", "dens"))
  
  iteration      <- x$iteration
  
  if (is.null(burnin)) {
    burnin       <- ceiling(iteration/2)
  }
  
  if (burnin + 2 >= iteration) {
    stop("burnin is too close to the number of MCMC steps performed.")
  }
  
  k              <- ncol(x$posterior)

  x              <- list(n.group     = x$n.group,
                         N           = x$N,
                         iteration   = x$iteration,
                         burnin      = burnin,
                         posterior   = x$posterior,
                         hyperparms  = x$hyperparms, 
                         accept.rate = x$accept.rate,
                         formula     = x$formula,
                         ctrl.mcmc   = x$ctrl.mcmc,
                         plot.type = plot.type)
  
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
  posterior[,k]   <- sqrt(posterior[,k])
  burnin          <- x$burnin
  iteration       <- x$iteration
  
  do.call("par", c(x[-(1:10)], list(...)))
  sim             <- x$plot.type == "sim"
  
  coln            <- colnames(posterior)
  coln[k - 1]     <- "Peer effect"
  if (sim) {
    for (i in 1:(k - 1)) {
      tmp         <- c(list(x = posterior[, i],
                            type = "l",
                            main = coln[i],
                            xlab = "",
                            ylab = "",
                            ...  = ...), x[-(1:11)])
      do.call("plot", tmp)
    }
    tmp           <- c(list(x = posterior[, k],
                             type = "l",
                          main = expression(sigma),
                          xlab = "",
                          ylab = "",
                          ...  = ...), x[-(1:11)])
    do.call("plot", tmp)
  } else {
    posterior     <- posterior[(burnin+1):iteration,,drop = FALSE]
    for (i in 1:(k - 1)) {
      tmp1        <- c(list(x = posterior[, i], ... = ...),  x[-(1:11)])

      tmp2        <- c(list(x = do.call("density", tmp1),
                            main = coln[i],
                            xlab = "",
                            ylab = "",
                            ...  = ...), x[-(1:11)])
      do.call("plot", tmp2)
    }
    tmp1          <- c(list(x = posterior[, k], ... = ...),  x[-(1:11)])
    
    tmp2          <- c(list(x = do.call("density", tmp1),
                          main = coln[i],
                          xlab = "",
                          ylab = "",
                          ...  = ...), x[-(1:11)])
    do.call("plot", tmp2)
  }
  

  par(mfrow = c(1, 1))

  invisible(x)
}

fsum      <- function(x, alpha) {
  c(mean(x), sd(x), quantile(x, probs = c(0.5*(1 - alpha), 0.5*(1 + alpha))))
} 


percent <- function(w, digits = 3, format = "f", ...) {
  paste0(formatC(100 * w, format = format, digits = digits, ...), "%")
}
