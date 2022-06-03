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

f.mc <- function(l, kappa){
  M       <- 20           # Number of groups
  N       <- rep(250,M)   # Group size
  Nsum    <- sum(N)
  Ncum    <- c(0, cumsum(N))
  
  # Parameters
  genzeta <- 1.5
  mu      <- -1.25
  sigma   <- 0.37
  K       <- 12
  P       <- 3
  
  distr   <- lapply(1:M, function(m){
    genz  <- rvMF(N[m], kappa*rvMF(1, rep(0, P)))
    gennu <- rnorm(N[m], mu, sigma)
    gend  <- N[m] * exp(gennu) * exp(mu + 0.5 * sigma ^ 2) * exp(logCpvMF(P, 0) - logCpvMF(P, genzeta))
    sim.dnetwork(nu = gennu, d = gend, zeta = genzeta, z = genz) 
  })
  
  X       <- cbind(X1 = rnorm(Nsum, 0, 2), X2 = rpois(Nsum, 2))
  theta   <- c(0.4, 2, 1, 1.5, 5, -3, 1)
  
  Ad      <- sim.network(dnetwork = distr, normalise = FALSE)
  G       <- norm.network(Ad)
  y       <- simsar(~ X, contextual = TRUE, Glist = G, theta = theta)
  Gy      <- y$Gy
  y       <- y$y
  
  # Compute W
  dG      <- sim.network(distr, normalise = TRUE)
  dGX     <- peer.avg(dG, X)
  dGdGX   <- peer.avg(dG, dGX)
  W       <- solve(crossprod(cbind(1, X, dGX, dGdGX))/Nsum)
  
  sest    <- smmSAR(y ~ X, dnetwork = distr, iv.power = 2L, W = W, smm.ctr  = list(print = F, R = 500),
                      fixed.effects = F, contextual = T)$estimates

  out     <- c(sest[-1], sest[1])
  cat(paste0(Sys.time(), " -- Iteration :", l), "\n")
  print(out)
  out
}

# Number of simulation
iteration     <- 500
kappa         <- c(0, 15, 30, 50)
n.kappa       <- length(kappa)

#######
set.seed(123)
out.mc        <- list()
for (x in 1:n.kappa) {
  # Construct cluster
  cl         <- makeCluster(25L)
  # After the function is run, close the cluster.
  on.exit(stopCluster(cl))
  # Register parallel backend
  registerDoParallel(cl)
  fn          <- paste0("log.kappa=", kappa[x], ".txt")
  if (file.exists(fn)) {file.remove(fn)}
  out.mc[[x]] <- foreach(l = 1:iteration, .combine = rbind, .packages = c("CDatanet", "PartialNetwork")) %dorng% 
    {sink(fn, append = TRUE); outx <- f.mc(l, kappa = kappa[x]); sink(); outx}
  save(out.mc, file = "mc.gx_unobserved.Rda")
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
  write.csv(results[[x]], file = paste0("~/Dropbox/Papers - In progress/Partial Network/Simulations/Monte Carlo/Results/Gx_unobserved_kappa=", kappa[x], ".csv"))
}

