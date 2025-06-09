##################################################################################
# This code replicates the Monte Carlo simulations when GX is observed and       #
# the network is estimated using nuclear ARD  (section 3.1).                     #
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
# parlambda is the tuning parameter for 
# nuclear estimation see Alidaee, H., E. Auerbach, and M. P. Leung (2020):
# “Recovering Network Structure from Aggregated Relational Data using Penalized Regression,” 
f.mc <- function(l, kappa){
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

  #loop over groups
  for (m in 1:M) {
    distr <- NULL
    while(is.null(distr)){try({
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
      
      #5 contruct ADR
      ARD               <- G %*% trait
      
      #Estimate the network distribution
      distr             <- accel_nuclear_gradient(inputs = t(trait), outputs = t(ARD), lambda = 600); ldistr[[m]] <- distr
    })}
    #True network row normalized
    W[[m]]   <- norm.network(G)
    
    #Covariates
    X[[m]]   <- cbind(rnorm(N[m],0,5),rpois(N[m],6))
    
    #True GX
    GX[[m]]  <- W[[m]] %*% X[[m]]
    
    #Y for section without contextual effects
    Y1[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1,N[m]), X[[m]])%*%beta + rnorm(N[m],0,se))
    GY1[[m]] <- W[[m]] %*% Y1[[m]]
    
    #Y for section with contextual effects
    Y2[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1, N[m]), X[[m]]) %*% beta +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY2[[m]] <- W[[m]] %*% Y2[[m]]
    
    #Compute instruments
    instr1 <- sim.IV(dnetwork = distr, X[[m]], Y1[[m]], replication = 1, power = 2)
    instr2 <- sim.IV(dnetwork = distr, X[[m]], Y2[[m]], replication = 1, power = 2)
    
    GY1c[[m]]   <- instr1[[1]]$G1y
    GX1c0[[m]]  <- instr1[[1]]$G1X[, , 1] 
    G2X1c0[[m]] <- instr1[[1]]$G1X[, , 2]  
    GX1c[[m]]   <- instr1[[1]]$G2X[, , 1]
    G2X1c[[m]]  <- instr1[[1]]$G2X[, , 2]
    
    
    GY2c[[m]]   <- instr2[[1]]$G1y
    GX2c0[[m]]  <- instr2[[1]]$G1X[, , 1] 
    G2X2c0[[m]] <- instr2[[1]]$G1X[, , 2]  
    GX2c[[m]]   <- instr2[[1]]$G2X[, , 1]
    G2X2c[[m]]  <- instr2[[1]]$G2X[, , 2]
  }
  
  # Concatenate M groups data
  # Y 
  Y1all     <- do.call("c", lapply(1:M, function(x) Y1[[x]]))
  Y2all     <- do.call("c", lapply(1:M, function(x) Y2[[x]]))
  
  # GY observed
  GY1all    <- do.call("c", lapply(1:M, function(x) GY1[[x]]))
  GY2all    <- do.call("c", lapply(1:M, function(x) GY2[[x]]))
  
  # GY constructed 
  GY1call   <- do.call("c", lapply(1:M, function(x) GY1c[[x]]))
  GY2call   <- do.call("c", lapply(1:M, function(x) GY2c[[x]]))
  
  # X
  Xall      <- do.call(rbind, lapply(1:M, function(x) X[[x]]))
  
  # GX observed
  GXall     <- do.call(rbind, lapply(1:M, function(x) GX[[x]]))
  
  # G^pX constructed
  # model without contextual effects
  # same draw
  GX1c0all   <- do.call(rbind, lapply(1:M, function(x) GX1c0[[x]]))
  G2X1c0all  <- do.call(rbind, lapply(1:M, function(x) G2X1c0[[x]]))
  # different draw
  GX1call    <- do.call(rbind, lapply(1:M, function(x) GX1c[[x]]))
  G2X1call   <- do.call(rbind, lapply(1:M, function(x) G2X1c[[x]]))
  
  # model with contextual effects
  # same draw
  GX2c0all   <- do.call(rbind, lapply(1:M, function(x) GX2c0[[x]]))
  G2X2c0all  <- do.call(rbind, lapply(1:M, function(x) G2X2c0[[x]]))
  # different draw
  GX2call    <- do.call(rbind, lapply(1:M, function(x) GX2c[[x]]))
  G2X2call   <- do.call(rbind, lapply(1:M, function(x) G2X2c[[x]]))
  
  ###### Estimation
  #Without contextual effects
  #if Gy is observed
  #iv
  sest1.1.1   <- summary(ivreg(Y1all ~ Xall + GY1all | Xall + GX1call + G2X1call), diagnostic=TRUE)
  lest1.1.1   <- c(sest1.1.1$coefficients[,1], sest1.1.1$diagnostics[,3])
  #smm
  W           <- solve(crossprod(cbind(1, Xall, GX1call, G2X1call))/sum(N))
  sest1.1.2   <- smmSAR(Y1all ~ Xall | GY1all, dnetwork = ldistr, iv.power = 2L, W = W, smm.ctr  = list(R = 1),
                        fixed.effects = F, contextual = F, cond.var = FALSE)$estimates
  lest1.1.2   <- c(sest1.1.2[-1], sest1.1.2[1])
  #if Gy is not observed
  #iv
  #Same draw
  sest1.2.1   <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1c0all + G2X1c0all), diagnostic=TRUE)
  lest1.2.1   <- c(sest1.2.1$coefficients[,1], sest1.2.1$diagnostics[,3], cor((GY1all-GY1call), cbind(GX1c0all,G2X1call)))
  #different draw
  sest1.2.2   <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1call + G2X1call), diagnostic=TRUE)
  lest1.2.2   <- c(sest1.2.2$coefficients[,1], sest1.2.2$diagnostics[,3], cor((GY1all-GY1call), cbind(GX1call,G2X1call)))
  
  #smm
  W           <- solve(crossprod(cbind(1, Xall, GX1call, G2X1call))/sum(N))
  sest1.2.3   <- smmSAR(Y1all ~ Xall, dnetwork = ldistr, iv.power = 2L, W = W, smm.ctr  = list(R = 1, S = 1, T = 1),
                        fixed.effects = F, contextual = F, cond.var = FALSE)$estimates
  lest1.2.3   <- c(sest1.2.3[-1], sest1.2.3[1])
  
  #With contextual effects
  #As GX is an explanatory variable, we distinguish whether GX is observed or not
  #GY is observed and GX is observed
  #iv
  sest2.1.1   <- summary(ivreg(Y2all ~ Xall + GXall + GY2all | Xall  + GXall + G2X2call), diagnostic = TRUE)
  lest2.1.1   <- c(sest2.1.1$coefficients[, 1], sest2.1.1$diagnostics[, 3])
  
  #smm
  W           <- solve(crossprod(cbind(1, Xall, GXall, peer.avg(sim.network(ldistr, normalise = TRUE), GXall)))/sum(N))
  sest2.1.2   <- smmSAR(Y2all ~ Xall | GY2all | GXall, dnetwork = ldistr, iv.power = 2L, W = W, smm.ctr  = list(R = 1),
                        fixed.effects = F, contextual = T, cond.var = FALSE)$estimates
  lest2.1.2   <- c(sest2.1.2[-1], sest2.1.2[1])
  
  #GY is not observed and GX is observed
  #iv
  sest2.2.1   <- summary(ivreg(Y2all ~ Xall + GXall + GX2c0all + GY2call | Xall  + GXall + GX2c0all + G2X2call), diagnostic = TRUE)
  lest2.2.1   <- c(sest2.2.1$coefficients[, 1], sest2.2.1$diagnostics[, 3])
  
  #smm
  W           <- solve(crossprod(cbind(1, Xall, GXall, peer.avg(sim.network(ldistr, normalise = TRUE), GXall)))/sum(N))
  sest2.2.2   <- smmSAR(Y2all ~ Xall || GXall, dnetwork = ldistr, iv.power = 2L, W = W, smm.ctr  = list(R = 1, S = 1, T = 1),
                        fixed.effects = F, contextual = T, cond.var = FALSE)$estimates
  lest2.2.2   <- c(sest2.2.2[-1], sest2.2.2[1])
  
  out         <-   c(lest1.1.1, lest1.1.2, lest1.2.1, lest1.2.2, lest1.2.3,
                    lest2.1.1, lest2.1.2, lest2.2.1, lest2.2.2)
  cat(paste0(Sys.time(), " -- Iteration :", l), "\n")
  print(out)
  out
}

# Number of simulation
iteration     <- 1000
kappa         <- c(0, 15, 30, 50)
n.kappa       <- length(kappa)

#######
set.seed(123)
out.mc        <- list()
for (x in 1:n.kappa) {
  # Construct cluster
  cl          <- makeCluster(10L)
  # After the function is run, close the cluster.
  on.exit(stopCluster(cl))
  # Register parallel backend
  registerDoParallel(cl)
  fn          <- paste0("log.kappa=", kappa[x], ".txt")
  if (file.exists(fn)) {file.remove(fn)}
  out.mc[[x]] <- foreach(l = 1:iteration, .combine = rbind, .packages = c("nuclearARD", "AER", "PartialNetwork")) %dorng% 
    {sink(fn, append = TRUE); outx <- f.mc(l, kappa = kappa[x]); sink(); outx}
  save(out.mc, file = "mc.gx_observed_Alidaee.Rda")
}

# the colnames
tmp1 <- c("Intercept", paste0("X",1:2), "alpha", "Weak", "Wu", "Sargan")
tmp2 <- c("Intercept", paste0("X",1:2), "alpha")
c1   <- paste0("No Con - GY obs - IV ", tmp1)
c2   <- paste0("No Con - GY obs - SMM ", tmp2)
c3   <- paste0("No Con - GY notobs - IV sam draw ", c(tmp1,"corGX1e","corGX2e","corGGX1e","corGGX2e"))
c4   <- paste0("No Con - GY notobs - IV dif draw ", c(tmp1,"corGX1e","corGX2e","corGGX1e","corGGX2e"))
c5   <- paste0("No Con - GY notobs - SMM ", tmp2)

tmp1 <- c("Intercept", paste0("X", 1:2), paste0("GX", 1:2), "alpha", "Weak", "Wu", "Sargan")
tmp2 <- c("Intercept", paste0("X", 1:2), paste0("GX", 1:2), "alpha")
c6   <- paste0("Wit Con - GY obs GX obs - IV ", tmp1)
c7   <- paste0("Wit Con - GY obs GX obs - SMM ", tmp2)
tmp1 <- c("Intercept", paste0("X", 1:2), paste0("GX", 1:2), paste0("GXc", 1:2), "alpha", "Weak", "Wu", "Sargan")
c8   <- paste0("Wit Con - GY notobs GX obs - IV ", tmp1)
c9   <- paste0("Wit Con - GY notobs GX obs - SMM ", tmp2)

# summary for all simulation using ARD
results        <-  lapply(1:n.kappa, function(x) {
  colnames(out.mc[[x]]) <- c(c1, c2, c3, c4, c5, c6, c7, c8, c9)
  t(apply(out.mc[[x]], 2, our.sum))})

print(results[[1]])
print(results[[2]])
print(results[[3]])
print(results[[4]])

for (x in 1:n.kappa) {
  write.csv(results[[x]], file = paste0("~/Dropbox/Papers - In progress/Partial Network/Simulations/Monte Carlo/Results/Gx_observed_Alidaee_kappa=", kappa[x], ".csv"))
}
