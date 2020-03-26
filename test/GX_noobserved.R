##################################################################################
# This code replicates the Monte Carlo simulations using Bayesian method.        #
# GX and Gy are not observed and the network distribution is used as             #
# prior information (section 4).                                                 #
#####################################Headline#####################################
rm(list = ls())
library(AER)                          # Install AER if not already done.
library(PartialNetwork)               # Install PartialNetwork if not already done.
library(doParallel)                   # To run the Monte Carlo in parallel
##################################################################################
# our summary function
our.sum <- function(x) {
  out <- c(mean(x, na.rm = TRUE),
           sd(x, na.rm = TRUE),
           quantile(x, 0.25, na.rm = TRUE),
           median(x, na.rm = TRUE),
           quantile(x, 0.75, na.rm = TRUE))
  names(out) <- c("Mean", "Sd.", "1st Qu.", "Median", "3rd Qu.")
  return(out)
}

# function to perform the simulation
# l stands for the l-th simulation

f.mc  <- function(l){
  cat("Iteration ", l, "\n")
  M          <- 100          # Number of groups
  N          <- rep(50,M)   # Group size
  lambda     <- 1
  
  beta       <- c(2,1,1.5)
  gamma      <- c(5,-3)
  alpha      <- 0.4
  se         <- 1
  
  ## generate data
  prior      <- list()
  # covariates
  X          <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
  # dependent variable
  y          <- c()
  
  #loop over groups
  for (m in 1:M) {
    Nm            <- N[m]
    #Generate link probabilities
    c             <- rnorm(Nm * Nm, 0, 1)
    
    priorm        <- matrix(exp(c / lambda) / (1 + exp(c / lambda)), Nm)
    diag(priorm)  <- 0
    prior[[m]]    <- priorm
    
    # true network
    Gm            <- matrix(runif(Nm^2), Nm, Nm) < priorm 
    # no self-link
    diag(Gm)      <- rep(0,Nm) 
    rsm           <- rowSums(Gm)
    rsm[rsm==0]   <- 1
    Gm            <- Gm/rsm 
    
    # rows index of group m
    r2            <- sum(N[1:m])
    r1            <- r2 - Nm + 1
    
    # contextual effects
    Xm            <- X[r1:r2,]
    GXm           <- Gm %*% Xm
    tmp           <- cbind(rep(1, Nm), Xm) %*% beta + GXm %*% gamma + rnorm(Nm, 0, se)
    y[r1:r2]      <- solve(diag(Nm) - alpha * Gm) %*% tmp
  }
  
  
  Kv           <- 2*ncol(X) + 1
  
  # prior information
  hyperparms   <- list("dnetwork" = prior,
                       "mutheta"  = rep(0,Kv),
                       "invsheta" = diag(Kv)/100,
                       "muzeta"   = 0,
                       "inszeta"  = 2,
                       "a"        = 4.2,
                       "b"        = (4.2 - 2)/se)
  ctrl         <- list("print.level" = 0)
  mcmcStep     <- 2000
  out          <- mcmcSAR(y ~ X | X, hyperparms = hyperparms, start = c(beta, gamma, alpha, se),
                          iteration = mcmcStep, ctrl.mcmc = ctrl)
  
  apply(out$posterior[500:mcmcStep,], 2, mean)
}


iteration      <- 1000
out.mc         <- mclapply(1:iteration, f.mc, mc.cores = 8L)


simu           <- t(do.call(cbind, out.mc))
colnames(simu) <- c("Intercept", paste0("X",1:2), paste0("GX",1:2), "alpha", "se2")

# results
results        <- t(apply(simu, 2, our.sum))
print(results)