#####################################Headline#####################################
rm(list = ls())
setwd("~/Dropbox/ARD/CppFiles")
library(ggplot2)                      # Install ggplot2 if not already done.
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
f.mc <- function(l){
  M         <- 100          # Number of groups
  N         <- rep(50,M)    # Group size
  # Parameters
  
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
  W    <- Y1   <- Y2   <- Y3  <- X    <- GY1 <- GY2  <- GY3  <- GX       <- list(M)
  GY1c <- GY2c <- GY3c <- GXc <- G2Xc <-GXc0 <-G2Xc0 <- G3Xc <- indexall <- list(M)
  
  #loop over group
  for (m in 1:M) {
    #Covariates
    X[[m]]   <- cbind(rnorm(N[m],0,5),rpois(N[m],6))
    
    #Fixed effects
    feffect  <- 0.3 * X[[m]][1, 1] + 0.3 * X[[m]][3, 2] - 1.8 * X[[m]][50, 2]
    
    # Network depending on X
    distr          <- matrix(0, N[m], N[m])
    for (i in 1:(N[m] - 1)) {
      for (j in (i + 1):N[m]) {
        distr[i,j] <- pnorm(-4.5 + abs(X[[m]][i,1] - X[[m]][j,1]) - 2*abs(X[[m]][i,2] - X[[m]][j,2]))
      }
    }
    distr          <- distr + t(distr)

    diag(distr)    <- 0
    
    #The complete graph
    G              <- sim.network(dnetwork = distr)
    
    #True network row normalized
    rs             <- rowSums(G)
    rs[rs == 0]    <- 1
    W[[m]]         <- G / rs
    
    #True GX
    GX[[m]]  <- W[[m]] %*% X[[m]]
    
    #Y for section without contextual effects
    Y1[[m]]  <- solve(diag(rep(N[m])) - alpha * W[[m]]) %*% (cbind(rep(1,N[m]), X[[m]])%*%beta + rnorm(N[m],0,se))
    GY1[[m]] <- W[[m]] %*% Y1[[m]]
    
    #Y for section with contextual effects
    Y2[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1, N[m]), X[[m]]) %*% beta +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY2[[m]] <- W[[m]] %*% Y2[[m]]
    
    #Y for section with Sub-Population Unobserved Fixed Effect
    Y3[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (feffect + X[[m]] %*% beta[-1] +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY3[[m]] <- W[[m]] %*% Y3[[m]]
    
    #Compute instruments
    instr1 <- sim.IV(dnetwork = distr, X[[m]], Y1[[m]], replication = 1, power = 2)
    instr2 <- sim.IV(dnetwork = distr, X[[m]], Y2[[m]], replication = 1, power = 2)
    instr3 <- sim.IV(dnetwork = distr, X[[m]], Y3[[m]], replication = 1, power = 2)
    
    GY1c[[m]]  <- instr1[[1]]$G1y
    GY2c[[m]]  <- instr2[[1]]$G1y
    GY3c[[m]]  <- instr3[[1]]$G1y
    
    GXc0[[m]]  <- instr1[[1]]$G1X[, , 1]
    G2Xc0[[m]] <- instr1[[1]]$G1X[, , 2]
    GXc[[m]]   <- instr1[[1]]$G2X[, , 1]
    G2Xc[[m]]  <- instr1[[1]]$G2X[, , 2]
  }
  
  # Concatenate M groups data
  
  # Y for section 3.1, 3.2 and 3.3
  Y1all     <- do.call("c", lapply(1:M, function(x) Y1[[x]]))
  Y2all     <- do.call("c", lapply(1:M, function(x) Y2[[x]]))
  Y3allm0   <- do.call("c", lapply(1:M, function(x) J %*% Y3[[x]]))
  
  
  # In section 3.1 or 3.3 GY may be observed
  GY1all    <- do.call("c", lapply(1:M, function(x) GY1[[x]]))
  GY2all    <- do.call("c", lapply(1:M, function(x) GY2[[x]]))
  GY3all    <- do.call("c", lapply(1:M, function(x) J %*% GY3[[x]]))
  
  # GY constructed for section 3.1
  GY1call   <- do.call("c", lapply(1:M, function(x) GY1c[[x]]))
  GY2call   <- do.call("c", lapply(1:M, function(x) GY2c[[x]]))
  GY3callm0 <- do.call("c", lapply(1:M, function(x) J %*% GY3c[[x]]))
  
  
  # Covariates
  Xall      <- do.call(rbind, lapply(1:M, function(x) X[[x]]))
  Xallm0    <- do.call(rbind, lapply(1:M, function(x) J %*% X[[x]]))
  
  # GX
  GXall     <- do.call(rbind, lapply(1:M, function(x) GX[[x]]))
  GXallm0   <- do.call(rbind, lapply(1:M, function(x) J %*% GX[[x]]))
  
  
  # G^pX constructed
  GXc0all   <- do.call(rbind, lapply(1:M, function(x) GXc0[[x]]))
  GXcall    <- do.call(rbind, lapply(1:M, function(x) GXc[[x]]))
  G2Xcall   <- do.call(rbind, lapply(1:M, function(x) G2Xc[[x]]))
  G2Xc0all  <- do.call(rbind, lapply(1:M, function(x) G2Xc0[[x]]))
  GXc0allm0 <- do.call(rbind, lapply(1:M, function(x) J %*% GXc0[[x]]))
  G2Xcallm0 <- do.call(rbind, lapply(1:M, function(x) J %*% G2Xc[[x]]))
  
  
  ###### Estimation
  #estimation section 3.1
  #if Gy is observed. We use GX constructed as instruments
  sest1.1.1   <- summary(ivreg(Y1all ~ Xall + GY1all | Xall +  GXcall), diagnostic=TRUE)
  lest1.1.1   <- c(sest1.1.1$coefficients[,1],sest1.1.1$diagnostics[,3])
  
  #if Gy is observed. We use GX and GGX constructed as instruments
  sest1.1.2   <- summary(ivreg(Y1all ~ Xall + GY1all | Xall +  GXcall + G2Xcall) , diagnostic=TRUE)
  lest1.1.2   <- c(sest1.1.2$coefficients[,1],sest1.1.2$diagnostics[,3])
  
  #if Gy is not observed. 
  #Same draw
  #We use GX constructed as instruments
  sest1.2.1.1 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GXc0all) , diagnostic=TRUE)
  lest1.2.1.1 <- c(sest1.2.1.1$coefficients[,1],sest1.2.1.1$diagnostics[,3], cor((GY1all-GY1call),GXc0all))
  
  #We use GX and GGX constructed as instruments
  sest1.2.1.2 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GXc0all + G2Xc0all), diagnostic=TRUE)
  lest1.2.1.2 <- c(sest1.2.1.2$coefficients[,1],sest1.2.1.2$diagnostics[,3], cor((GY1all-GY1call),cbind(GXc0all,G2Xcall)))
  
  #different draw
  #We use GX constructed as instruments
  sest1.2.2.1 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GXcall), diagnostic=TRUE)
  lest1.2.2.1 <- c(sest1.2.2.1$coefficients[,1],sest1.2.2.1$diagnostics[,3], cor((GY1all-GY1call),GXcall))
  
  #We use GX and GGX constructed as instruments
  sest1.2.2.2 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GXcall + G2Xcall), diagnostic=TRUE)
  lest1.2.2.2 <- c(sest1.2.2.2$coefficients[,1],sest1.2.2.2$diagnostics[,3], cor((GY1all-GY1call),cbind(GXcall,G2Xcall)))
  
  #estimation section 3.2
  #GY is observed
  sest2.1     <- summary(ivreg(Y2all ~ Xall + GXall + GY2all | Xall  + GXall + G2Xcall), diagnostic = TRUE)
  lest2.1     <- c(sest2.1$coefficients[, 1], sest2.1$diagnostics[, 3])
  
  #GY is not observed
  sest2.2     <- summary(ivreg(Y2all ~ Xall + GXall +  GXc0all + GY2call | Xall  + GXall + GXc0all + G2Xcall), diagnostic = TRUE)
  lest2.2     <- c(sest2.2$coefficients[, 1], sest2.2$diagnostics[, 3])
  
  # estimation section 3.3
  sest3.1     <- summary(ivreg(Y3allm0 ~ Xallm0 + GXallm0 + GY3all | Xallm0 + GXallm0 + G2Xcallm0), diagnostic = TRUE)
  lest3.1     <- c(sest3.1$coefficients[, 1], sest3.1$diagnostics[, 3])
  
  #Y is not observed
  sest3.2     <- summary(ivreg(Y3allm0 ~ Xallm0 + GXallm0 + GXc0allm0 + GY3callm0 | Xallm0 + GXallm0 + GXc0allm0 + G2Xcallm0), diagnostic = TRUE)
  lest3.2     <- c(sest3.2$coefficients[, 1], sest3.2$diagnostics[, 3])
  
  cat("Iteration ", l, "\n")
  c(lest1.1.1, lest1.1.2, lest1.2.1.1, lest1.2.1.2, lest1.2.2.1,
    lest1.2.2.2, lest2.1, lest2.2, lest3.1, lest3.2)
}

set.seed(123)

# Number of simulation
iteration <- 1000
out.mc        <- mclapply(1:iteration, f.mc, mc.cores = 8L)

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
print(results)