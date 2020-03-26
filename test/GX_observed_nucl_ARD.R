##################################################################################
# This code replicates the Monte Carlo simulations when GX is observed and       #
# the network is estimated using nuclear ARD.                                    #
#####################################Headline#####################################
rm(list = ls())
library(AER)                          # Install AER if not already done.
library(PartialNetwork)               # Install PartialNetwork if not already done.
library(nuclearARD)                   # Install nuclearARD if not already done.
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
# parlambda is the tuning parameter for 
# nuclear estimation see Alidaee, H., E. Auerbach, and M. P. Leung (2020):
# “Recovering Network Structure from Aggregated Relational Data using Penalized Regression,” 
fsim <- function(l, parlambda){
  M          <- 20           # Number of groups
  N          <- rep(250,M)   # Group size
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
  
  P1        <- matrix(1, N[1], N[1])
  diag(P1)  <- 0
  rs        <- rowSums(P1)
  P1        <- P1 / rs
  J         <- diag(N[1]) - P1
  
  
  # Some useful variables
  W     <- Y1    <- Y2   <- Y3    <- X     <- GY1   <- GY2    <- GY3      <- GX  <- list(M)
  GY1c  <- GY2c  <- GY3c <- GX1c  <- G2X1c <-GX1c0  <- G2X1c0 <- indexall <- list(M)
  GX2c  <- G2X2c <-GX2c0 <-G2X2c0 <- GX3c  <- G2X3c <- GX3c0  <-G2X3c0    <- list(M)
  
  #loop over groups
  for (m in 1:M) {
    print(paste("Iteration :", l," -- Group:", m))
    #1- Generate z
    genz  <- rvMF(N[m], rep(0, P))
    #2- Genetate nu  from a Normal distribution with parameters mu and sigma
    gennu <- rnorm(N[m], mu, sigma)
    #3- Generate a graph G
    #Before, lets's compute d
    gend  <- N[m] * exp(gennu) * exp(mu + 0.5 * sigma ^ 2) * exp(logCpvMF(P, 0) - logCpvMF(P, genzeta))
    
    #Link probabilities
    Probabilities     <- sim.dnetwork(nu = gennu, d = gend, zeta = genzeta, z = genz) 
    
    #The complete graph
    G                 <- sim.network(dnetwork = Probabilities)
    
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
    distr    <- accel_nuclear_gradient(inputs = t(trait), outputs = t(ARD), lambda = parlambda)
    
    #True network row normalized
    rs               <- rowSums(G)
    rs[rs == 0]      <- 1
    W[[m]]           <- G / rs
    
    #Covariates
    X[[m]]   <- cbind(rnorm(N[m],0,5),rpois(N[m],6))
    
    #Fixed effects
    feffect  <- 0.3 * X[[m]][1, 1] + 0.3 * X[[m]][3, 2] - 1.8 * X[[m]][50, 2]
    
    #True GX
    GX[[m]]  <- W[[m]] %*% X[[m]]
    
    #Y for section without contextual effects
    Y1[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1,N[m]), X[[m]])%*%beta + rnorm(N[m],0,se))
    GY1[[m]] <- W[[m]] %*% Y1[[m]]
    
    #Y for section with contextual effects
    Y2[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1, N[m]), X[[m]]) %*% beta +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY2[[m]] <- W[[m]] %*% Y2[[m]]
    
    #Y for section with Sub-Population Unobserved Fixed Effect
    Y3[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (feffect + X[[m]] %*% beta[-1] +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY3[[m]] <- W[[m]] %*% Y3[[m]]
    
    
    #Estimate the network distribution
    distr  <- sim.network(distr)
    
    #Compute instruments
    instr1 <- sim.IV(dnetwork = distr, X[[m]], Y1[[m]], replication = 1, power = 2)
    instr2 <- sim.IV(dnetwork = distr, X[[m]], Y2[[m]], replication = 1, power = 2)
    instr3 <- sim.IV(dnetwork = distr, X[[m]], Y3[[m]], replication = 1, power = 2)
    
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
    
    
    GY3c[[m]]   <- instr3[[1]]$G1y
    GX3c0[[m]]  <- instr3[[1]]$G1X[, , 1] 
    G2X3c0[[m]] <- instr3[[1]]$G1X[, , 2]  
    GX3c[[m]]   <- instr3[[1]]$G2X[, , 1]
    G2X3c[[m]]  <- instr3[[1]]$G2X[, , 2]
  }
  
  # Concatenate M groups data
  # Y 
  Y1all     <- do.call("c", lapply(1:M, function(x) Y1[[x]]))
  Y2all     <- do.call("c", lapply(1:M, function(x) Y2[[x]]))
  Y3all     <- do.call("c", lapply(1:M, function(x) J %*% Y3[[x]]))
  
  
  # GY observed
  GY1all    <- do.call("c", lapply(1:M, function(x) GY1[[x]]))
  GY2all    <- do.call("c", lapply(1:M, function(x) GY2[[x]]))
  GY3all    <- do.call("c", lapply(1:M, function(x) J %*% GY3[[x]]))
  
  # GY constructed 
  GY1call   <- do.call("c", lapply(1:M, function(x) GY1c[[x]]))
  GY2call   <- do.call("c", lapply(1:M, function(x) GY2c[[x]]))
  GY3call   <- do.call("c", lapply(1:M, function(x) J %*% GY3c[[x]]))
  
  
  # X
  Xall      <- do.call(rbind, lapply(1:M, function(x) X[[x]]))
  X3all     <- do.call(rbind, lapply(1:M, function(x) J %*% X[[x]]))
  
  # GX observed
  GXall     <- do.call(rbind, lapply(1:M, function(x) GX[[x]]))
  GX3all    <- do.call(rbind, lapply(1:M, function(x) J %*% GX[[x]]))
  
  
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
  
  # model with fixed effects
  # same draw
  GX3c0all   <- do.call(rbind, lapply(1:M, function(x) J %*% GX3c0[[x]]))
  G2X3c0all  <- do.call(rbind, lapply(1:M, function(x) J %*% G2X3c0[[x]]))
  # different draw
  GX3call    <- do.call(rbind, lapply(1:M, function(x) J %*% GX3c[[x]]))
  G2X3call   <- do.call(rbind, lapply(1:M, function(x) J %*% G2X3c[[x]]))
  
  
  ###### Estimation
  #estimation section 3.1
  #if Gy is observed. We use GX constructed as instruments
  sest1.1.1   <- summary(ivreg(Y1all ~ Xall + GY1all | Xall +  GX1call), diagnostic=TRUE)
  lest1.1.1   <- c(sest1.1.1$coefficients[,1],sest1.1.1$diagnostics[,3])
  
  #if Gy is observed. We use GX and GGX constructed as instruments
  sest1.1.2   <- summary(ivreg(Y1all ~ Xall + GY1all | Xall +  GX1call + G2X1call) , diagnostic=TRUE)
  lest1.1.2   <- c(sest1.1.2$coefficients[,1],sest1.1.2$diagnostics[,3])
  
  #if Gy is not observed. 
  #Same draw
  #We use GX constructed as instruments
  sest1.2.1.1 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1c0all) , diagnostic=TRUE)
  lest1.2.1.1 <- c(sest1.2.1.1$coefficients[,1],sest1.2.1.1$diagnostics[,3])
  
  #We use GX and GGX constructed as instruments
  sest1.2.1.2 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1c0all + G2X1c0all), diagnostic=TRUE)
  lest1.2.1.2 <- c(sest1.2.1.2$coefficients[,1],sest1.2.1.2$diagnostics[,3])
  
  #different draw
  #We use GX constructed as instruments
  sest1.2.2.1 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1call), diagnostic=TRUE)
  lest1.2.2.1 <- c(sest1.2.2.1$coefficients[,1],sest1.2.2.1$diagnostics[,3])
  
  #We use GX and GGX constructed as instruments
  sest1.2.2.2 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1call + G2X1call), diagnostic=TRUE)
  lest1.2.2.2 <- c(sest1.2.2.2$coefficients[,1],sest1.2.2.2$diagnostics[,3])
  
  #estimation section 3.2
  #GY is observed
  sest2.1     <- summary(ivreg(Y2all ~ Xall + GXall + GY2all | Xall  + GXall + G2X2call), diagnostic = TRUE)
  lest2.1     <- c(sest2.1$coefficients[, 1], sest2.1$diagnostics[, 3])
  
  #GY is not observed
  sest2.2     <- summary(ivreg(Y2all ~ Xall + GXall +  GX2c0all + GY2call | Xall  + GXall + GX2c0all + G2X2call), diagnostic = TRUE)
  lest2.2     <- c(sest2.2$coefficients[, 1], sest2.2$diagnostics[, 3])
  
  # estimation section 3.3
  #Gy is observed
  sest3.1     <- summary(ivreg(Y3all ~ X3all  + GX3all + GY3all | X3all + GX3all + G2X3call), diagnostic = TRUE)
  lest3.1     <- c(sest3.1$coefficients[, 1], sest3.1$diagnostics[, 3])
  
  #Gy is not observed
  sest3.2     <- summary(ivreg(Y3all ~ X3all + GX3all + GX3c0all + GY3call | X3all + GX3all + GX3c0all + G2X3call), diagnostic = TRUE)
  lest3.2     <- c(sest3.2$coefficients[, 1], sest3.2$diagnostics[, 3])
  
  c(lest1.1.1, lest1.1.2, lest1.2.1.1, lest1.2.1.2, lest1.2.2.1,
    lest1.2.2.2, lest2.1, lest2.2, lest3.1, lest3.2)
}

# monte carlo function for each parlambda
f.mc <- function(iteration, parlambda) {
  out.mc        <- mclapply(1:iteration, function(w) fsim(w, parlambda), mc.cores = 8L)
  # simu as m matrix
  simu          <- t(do.call(cbind, out.mc))
  
  # the colnames
  tmp <- c("Intercept",paste0("X",1:2),"alpha","Weak","Wu","Sargan")
  c1  <- paste0("No Con - GY obs - ins GX ", tmp)
  c2  <- paste0("No Con - GY obs - ins GX GGX ", tmp)
  c3  <- paste0("No Con - GY notobs - ins GX - sam draw ", c(tmp,"corGX1e","corGX2e"))
  c4  <- paste0("No Con - GY notobs - ins GX GGX - sam draw ", c(tmp,"corGX1e","corGX2e","corGGX1e","corGGX2e"))
  c5  <- paste0("No Con - GY notobs - ins GX GGX - dif draw ", c(tmp,"corGX1e","corGX2e"))
  c6  <- paste0("No Con - GY notobs - ins GX GGX - dif draw ", c(tmp,"corGX1e","corGX2e","corGGX1e","corGGX2e"))
  
  tmp <- c("Intercept", paste0("X", 1:2), paste0("GX", 1:2), "alpha", "Weak", "Wu", "Sargan")
  c7  <- paste0("Wit Con - GY obs ", tmp)
  c9  <- paste0("Fix eff - GY obs ", tmp)
  tmp <- c("Intercept", paste0("X", 1:2), paste0("GX", 1:2), paste0("GXc", 1:2), "alpha", "Weak", "Wu", "Sargan")
  c8  <- paste0("Wit Con - GY no obs ", tmp)
  c10 <- paste0("Fix eff - GY obs ", tmp)
  
  colnames(simu) <- c(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
  
  # summary for all simulation using ARD
  results        <- t(apply(simu, 2, our.sum))
  results
}

set.seed(123)

# Number of simulation
iteration <- 1000
out1      <- f.mc(iteration, 200)
out2      <- f.mc(iteration, 600)
out3      <- f.mc(iteration, 1374)

print(out1)
print(out2)
print(out3)