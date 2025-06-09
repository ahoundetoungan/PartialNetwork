##################################################################################
# This code replicates the Monte Carlo simulations when GX is not observed and   #
# the network is estimated using latent  space model with ARD (section 3.1).     #
#####################################Headline#####################################
rm(list = ls())
library(doParallel)                   # To run the Monte Carlo in parallel
library(foreach)                      # To run the Monte Carlo in parallel
library(doRNG)                        # To run the Monte Carlo in parallel
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
# kappa is the concentration parameter
f.mc  <- function(l, kappa){
  M         <- 20           # Number of groups
  N         <- rep(250,M)   # Group size
  # Parameters
  genzeta   <- 1.5
  mu        <- -1.25
  sigma     <- 0.37
  K         <- 12
  P         <- 3
  
  beta      <- c(2,1,1.5)
  gamma     <- c(5,-3)
  alpha     <- 0.4
  se        <- 1
  
  
  # Some useful variables
  W     <- Y1    <- Y2     <- X        <- GY1   <- GY2   <- GX <-  GY1c   <- GY2c   <- GX1c   <- list(M)
  G2X1c <- GX1c0 <- G2X1c0 <- indexall <- GX2c  <- G2X2c <-GX2c0 <- G2X2c0 <- ldistr <- list(M)
  
  #loop over group to estimate dnetwork
  for (m in 1:M) {
    #1- Generate z
    genz  <- rvMF(N[m], kappa*rvMF(1, rep(0, P)))
    #2- Genetate nu  from a Normal distribution with parameters mu and sigma
    gennu <- rnorm(N[m], mu, sigma)
    #3- Generate a graph G
    #Before, let's compute d
    gend  <- N[m] * exp(gennu) * exp(mu + 0.5 * sigma ^ 2) * exp(logCpvMF(P, 0) - logCpvMF(P, genzeta))
    
    #Link probabilities
    Probabilities     <- sim.dnetwork(nu = gennu, d = gend, zeta = genzeta, z = genz) 
    
    #The complete graph
    G                 <- sim.network(Probabilities)
    
    #4a Generate vk
    genv              <- rvMF(K, rep(0, P))
    
    #fix some vk distant
    genv[1, ]         <- c(1, 0, 0)
    genv[2, ]         <- c(0, 1, 0)
    genv[3, ]         <- c(0, 0, 1)
    
    #4b set eta
    geneta            <- abs(rnorm(K, 4, 1))
    
    #4c Build Features matrix
    densityatz        <- matrix(0, N[m], K)
    for (k in 1:K) {
      densityatz[, k] <- dvMF(genz, genv[k, ] * geneta[k])
    }
    
    trait             <- matrix(0, N[m], K)
    
    NK                <- floor(runif(K, 0.8, 0.95) * colSums(densityatz) / unlist(lapply(1:K, function(w){max(densityatz[,w])}))) 
    
    for (k in 1:K) {
      trait[,k]       <- rbinom(N[m], 1, NK[k] * densityatz[,k] / sum(densityatz[,k]))
    } 
    
    #5 construct ADR
    ARD               <- G %*% trait
    
    #generate b
    genb              <- numeric(K)
    for (k in 1:K) {
      genb[k]         <- sum(G[, trait[, k] == 1]) / sum(G) + 1e-8
    }
    ###### Fit Breza et al. (2017) parameter #####
    # We should fix one bk
    vfixcolumn       <- 1:5
    bfixcolumn       <- c(1, 3, 5, 7, 11)
    
    # Initialization
    # We initialize using the true parameter in order to run just few MCMC steps
    start      <-list("z"    = genz,
                      "v"    = genv,
                      "d"    = gend,
                      "b"    = genb,
                      "eta"  = geneta,
                      "zeta" = genzeta)
    
    # fit the parameters
    # we consider the model with zeta = 1
    Gparms   <- mcmcARD(ARD, trait, start, vfixcolumn, bfixcolumn, iteration = 1000L,
                        ctrl.mcmc = list(print = FALSE))
    
    #Estimate the network distribution
    distr    <- fit.dnetwork(Gparms, print = FALSE)$dnetwork; ldistr[[m]] <- distr
    
    #True network row normalized
    W[[m]]   <- norm.network(G)
    
    #Covariates
    X[[m]]   <- cbind(rnorm(N[m],0,5),rpois(N[m],6))
    
    #True GX
    GX[[m]]  <- W[[m]] %*% X[[m]]
    
    #Y for section with contextual effects
    Y2[[m]]  <- solve(diag(N[m]) - alpha * W[[m]], cbind(rep(1, N[m]), X[[m]]) %*% beta +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY2[[m]] <- W[[m]] %*% Y2[[m]]
  }
  
  # Concatenate M groups data
  # Y 
  Y2all     <- do.call("c", lapply(1:M, function(x) Y2[[x]]))
  
  # GY observed
  GY2all    <- do.call("c", lapply(1:M, function(x) GY2[[x]]))
  
  # X
  Xall      <- do.call(rbind, lapply(1:M, function(x) X[[x]]))
  
  # Compute W
  dG        <- sim.network(ldistr, normalise = TRUE)
  dGXall    <- peer.avg(dG, Xall)
  dGdGXall  <- peer.avg(dG, dGXall)
  W         <- solve(crossprod(cbind(1, Xall, dGXall, dGdGXall))/sum(N))
  
  # #GY is observed and GX is not observed
  # #smm
  # sest2.3   <- smmSAR(Y2all ~ Xall | GY2all, dnetwork = ldistr, iv.power = 2L, W = W, smm.ctr  = list(R = 15, S = 15),
  #                     fixed.effects = F, contextual = T)$estimates
  # lest2.3   <- c(sest2.3[-1], sest2.3[1])

  #GY is not observed and GX is not observed
  #smm
  sest2.4   <- smmSAR(Y2all ~ Xall, dnetwork = ldistr, iv.power = 2L, W = W, smm.ctr  = list(R = 500),
                        fixed.effects = F, contextual = T)$estimates
  out       <- c(sest2.4[-1], sest2.4[1])
  cat(paste0(Sys.time(), " -- Iteration :", l), "\n")
  print(out)
  out
}

# Number of simulation
iteration     <- 500
kappa         <- c(15, 0, 30, 50)
n.kappa       <- length(kappa)

#######
set.seed(123)
out.mc        <- list()
for (x in 1:n.kappa) {
  # Construct cluster
  cl         <- makeCluster(5L)
  # After the function is run, close the cluster.
  on.exit(stopCluster(cl))
  # Register parallel backend
  registerDoParallel(cl)
  fn          <- paste0("log.kappa=", kappa[x], ".txt")
  if (file.exists(fn)) {file.remove(fn)}
  out.mc[[x]] <- foreach(l = 1:iteration, .combine = rbind, .packages = c("AER", "PartialNetwork")) %dorng% 
    {sink(fn, append = TRUE); outx <- f.mc(l, kappa = kappa[x]); sink(); outx}
  save(out.mc, file = "mc.gx_unobserved_Breza.Rda")
}

# the colnames
c10  <- paste0("Wit Con - GY notobs GX notobs - SMM ", c("Intercept", paste0("X", 1:2), paste0("GX", 1:2), "alpha"))

# summary for all simulation using ARD
results        <-  lapply(1:n.kappa, function(x) {
  colnames(out.mc[[x]]) <- c10
  t(apply(out.mc[[x]], 2, our.sum))})

print(results[[1]])
print(results[[2]])
print(results[[3]])
print(results[[4]])

for (x in 1:n.kappa) {
  write.csv(results[[x]], file = paste0("~/Dropbox/Papers - In progress/Partial Network/Simulations/Monte Carlo/Results/Gx_unobserved_Breza_kappa=", kappa[x], ".csv"))
}

