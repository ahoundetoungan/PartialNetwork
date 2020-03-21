#####################################Headline#####################################
## Required package Rcpp RcppArmadillo BH
rm(list = ls())
#setwd("~/Dropbox/ARD/Simulations/Monte Carlo") # for my laptop
#setwd("~/ARD/Simulations/Monte Carlo") # for the cluster
library(ggplot2)
library(AER)
library(nuclearARD)
library(doParallel)
library(PartialNetwork)

library(Rmpi)
##################################################################################

####################################Simulations###################################
set.seed(123)


# CpvMF 
CpvMF   <- function(p, k) {
  exp(logCpvMF(p, k))
}

#1000 MC
f <- function(l, parlambda){
  M <- 20            # Number of groups
  N <- rep(250, M)   # Group size
  
  
  # Parameters
  genzeta <- 1.5
  mu      <- -1.25
  sigma   <- 0.37
  K       <- 12
  P       <- 3
  
  beta  <- c(2,1,1.5)
  gamma <- c(5,-3)
  alpha <- 0.4
  se    <- 1
  
  P1       <- matrix(1, N[1], N[1])
  diag(P1) <- 0
  rs       <- rowSums(P1)
  P1       <- P1 / rs
  J        <- diag(N[1]) - P1
  
  # Some useful variables
  W <- Y1 <- Y2 <- Y3 <- X <- GY1 <- GY2 <- GY3 <- GX <-
    GY1c <- GY2c <- GY3c <- GXc <- G2Xc <-GXc0<-G2Xc0 <- G3Xc <- indexall <- list(M)
  
  #loop for group
  #We can run this loop in parallel to make the code faster
  for (m in 1:M) {
    #1- Generate z
    genz <- rvMF(N[m], rep(0, P))
    #2- Genetate nu  from a Normal distribution with parameters mu and sigma
    gennu <- rnorm(N[m], mu, sigma)
    #3- Generate a graph G
    #Before, lets's compute d
    gend <- N[m] * exp(gennu) * exp(mu + 0.5 * sigma ^ 2) * (CpvMF(P, 0) / CpvMF(P, genzeta))
    #Link probabilities
    #The function Prob also prints the maximum of link probabilities. If éthe maximum is greater than 1, please handle
    #the parameter genzeta, mu  or sigma. You may decrease mu to -1.37 or -1.4; or set genzeta to 0.9 for example.
    #Moreover, choose parameters so that the  maximum is not too low.
    
    Probabilities <- Prob(gennu, gend, genzeta, genz) 
    
    #The complete graph
    #G<-matrix(rep(0,N*N),N)
    G <- Graph(Probabilities)
    #probtograph(Probabilities,G,gend,N)
    
    #4a Generate vk
    genv <- rvMF(K, rep(0, P))
    #fix some vk distant
    genv[1, ] = c(1, 0, 0)
    genv[2, ] = c(0, 1, 0)
    genv[3, ] = c(0, 0, 1)
    
    
    #4b set eta
    geneta <- abs(rnorm(K, 4, 1))
    #4c Build Features matrix
    
    #compure density for each z
    densityatz <- matrix(0, N[m], K)
    for (k in 1:K) {
      densityatz[, k] <- dvMF(genz, genv[k, ] * geneta[k])
    }
    
    trait <- matrix(0, N[m], K)
    
    #The method to build features is not still defined. The are making some tests
    # Choose only one pattern to test the code
    
    #Pattern 1
    #Individual i has a feature k with the probability  f(zi,vk)/max(f(zi,v1),f(zi,v2),...,f(zi,vK))
    
    #for(i in 1:N){
    #  rowmax=max(densityatz[i,])
    #  for(k in 1:K){
    #    trait[i,k]=(runif(1)<densityatz[i,k]/rowmax)
    #  }
    #}
    
    #Pattern 2
    #For each k, take the Nk individuals closest to vk (that is, Nk highest f(zi,vk))
    #Nk is sampled from the uniform distribution on [0.04*N,0.12*N]. Theses values was calibrated so that the parameter converge
    #towards to simulated values. 5% to 12% of individuals have a features. This rate may seem low, but remember that the
    #features are rare trait so that the Bernoulli approximation to Poisson makes sense.
    
    #For now, this pattern seem to be the best among those tested
    #for (k in 1:K) {
    #  trait[, k] = (densityatz[, k] > sort(densityatz[, k], decreasing = T)[runif(1, 0.05 *
    #                                                                                N[m], 0.25 * N[m])])
    #}
    
    
    #Pattern 3
    #Each i has a feature k with the probability Nk*f(zi,vk)/sum_{k=1,2,...,K}(f(zi,vk))
    #Nk is sampled from the uniform distribution on [0.04*N,0.12*N].
    #This pattern is close of the patter 2. It is the random case the pattern 2
    #It also offers better results
    
    NK <- floor(runif(K,0.8,0.95)*colSums(densityatz)/unlist(lapply(1:K, function(w){max(densityatz[,w])}))) # added VB
    
    for (k in 1:K) { # added VB
      trait[,k]<-rbinom(N[m],1,NK[k]*densityatz[,k]/sum(densityatz[,k])) # added VB
    } # added VB
    
    
    #print a percentage of peaople having a trait
    colSums(trait) * 100 / N[m]
    
    #5 contruct ADR
    ARD <- G %*% trait
 
    #Estimate the network distribution
    
    distr  <- accel_nuclear_gradient(inputs = t(trait), outputs = t(ARD), lambda = parlambda)
    
    #True network row normalized
    W[[m]] <- G / rowSums(G)
    
    #Covariates
    X[[m]] <- cbind(rnorm(N[m],0,5),rpois(N[m],6))
    
    #Fixed effects
    feffect <- 0.3 * X[[m]][1, 1] + 0.3 * X[[m]][3, 2] - 1.8 * X[[m]][50, 2]
    
    #True GX
    GX[[m]] <- W[[m]] %*% X[[m]]
    
    #Y for sections 3.1 without contextual effects
    Y1[[m]] <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1,N[m]), X[[m]])%*%beta + rnorm(N[m],0,se))
    GY1[[m]] <- W[[m]] %*% Y1[[m]]
    
    #Y for section 3.2 with contextual effects
    Y2[[m]] <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1, N[m]), X[[m]]) %*% beta +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY2[[m]] <- W[[m]] %*% Y2[[m]]
    
    #Y for section 3.3 Sub-Population Unobserved Fixed Effect
    Y3[[m]] <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (feffect + X[[m]] %*% beta[-1] +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY3[[m]] <- W[[m]] %*% Y3[[m]]
    
    

    #Compute instruments
    # The function has 5 arguments
    # The network distribution distr and the explanatory variables X are required
    # y, S and pow are optional
    # S (by default = 2) is the number of draws of G needed for GX such that the first draw corresponds to the draw G used to compute GY
    # pow is used to compute GX, G^2X, ..., G^(pow)X
    
    # The output of this function is as follow
    # $GY  estimates GY by drawin one G from the network distribution and multiplies it by Y
    # $GX is a list of S instrument matrices, Each index in the list corresponds to one draw of G to compute GX. The first draw is also used to
    # compute GY
    # Each matrix instrument is cube (three dimensions) (N,K,pow)
    # the first axe indexes the individual, the second axe index the variable,
    # and the fird axe the power (GX, G^2X, G^3X)
    # For example $GX[[2]][,,3] is G^3X where G is the second draw of G
    
    instr1 <- instruments(distr, X[[m]], Y1[[m]], pow = 2)
    instr2 <- instruments(distr, X[[m]], Y2[[m]], pow = 2)
    instr3 <- instruments(distr, X[[m]], Y3[[m]], pow = 2)
    
    GY1c[[m]] <- instr1$GY
    GY2c[[m]] <- instr2$GY
    GY3c[[m]] <- instr3$GY
    
    GXc0[[m]] <- instr1$GX[[1]][, , 1]
    G2Xc0[[m]] <- instr1$GX[[1]][, , 2]
    GXc[[m]] <- instr1$GX[[2]][, , 1]
    G2Xc[[m]] <- instr1$GX[[2]][, , 2]
    print(paste(l, m))
  }
  
  # Concatenate M groups data
  
  # Y for section 3.1, 3.2 and 3.3
  Y1all <- do.call("c", lapply(1:M, function(x) {
    Y1[[x]]
  }))
  Y2all <- do.call("c", lapply(1:M, function(x) {
    Y2[[x]]
  }))
  Y3allm0 <- do.call("c", lapply(1:M, function(x) {
    J %*% Y3[[x]]
  }))
  
  
  # In section 3.1 or 3.3 GY may be observed
  GY1all <- do.call("c", lapply(1:M, function(x) {
    GY1[[x]]
  }))
  GY2all <- do.call("c", lapply(1:M, function(x) {
    GY2[[x]]
  }))
  GY3all <- do.call("c", lapply(1:M, function(x) {
    J %*% GY3[[x]]
  }))
  
  # GY constructed for section 3.1
  GY1call <- do.call("c", lapply(1:M, function(x) {
    GY1c[[x]]
  }))
  GY2call <- do.call("c", lapply(1:M, function(x) {
    GY2c[[x]]
  }))
  GY3callm0 <- do.call("c", lapply(1:M, function(x) {
    J %*% GY3c[[x]]
  }))
  
  
  # Covariates
  Xall <- do.call(rbind, lapply(1:M, function(x) {
    X[[x]]
  }))
  Xallm0 <- do.call(rbind, lapply(1:M, function(x) {
    J %*% X[[x]]
  }))
  
  # GX
  GXall <- do.call(rbind, lapply(1:M, function(x) {
    GX[[x]]
  }))
  GXallm0 <- do.call(rbind, lapply(1:M, function(x) {
    J %*% GX[[x]]
  }))
  
  
  # G^pX constructed
  GXc0all <- do.call(rbind, lapply(1:M, function(x) {
    GXc0[[x]]
  }))
  GXcall <- do.call(rbind, lapply(1:M, function(x) {
    GXc[[x]]
  }))
  G2Xcall <- do.call(rbind, lapply(1:M, function(x) {
    G2Xc[[x]]
  }))
  G2Xc0all <- do.call(rbind, lapply(1:M, function(x) {
    G2Xc0[[x]]
  }))
  GXc0allm0 <- do.call(rbind, lapply(1:M, function(x) {
    J %*% GXc0[[x]]
  }))
  G2Xcallm0 <- do.call(rbind, lapply(1:M, function(x) {
    J %*% G2Xc[[x]]
  }))
  
  
  ###### Estimation
  # estimation section 3.1
  # if Gy is observed. We use GX constructed as instruments
  est1.1.1<-ivreg(Y1all ~ Xall + GY1all | Xall +  GXcall)
  sest1.1.1<-summary(est1.1.1, diagnostic=TRUE)
  lest1.1.1 <- c(sest1.1.1$coefficients[,1],sest1.1.1$diagnostics[,3])
  
  # if Gy is observed. We use GX and GGX constructed as instruments
  est1.1.2<-ivreg(Y1all ~ Xall + GY1all | Xall +  GXcall + G2Xcall)
  sest1.1.2<-summary(est1.1.2, diagnostic=TRUE)
  lest1.1.2 <- c(sest1.1.2$coefficients[,1],sest1.1.2$diagnostics[,3])

  
  # if Gy is not observed. 
  #Same draw
  #We use GX constructed as instruments
  est1.2.1.1<-ivreg(Y1all ~ Xall + GY1call | Xall +  GXc0all)
  sest1.2.1.1<-summary(est1.2.1.1, diagnostic=TRUE)
  lest1.2.1.1 <- c(sest1.2.1.1$coefficients[,1],sest1.2.1.1$diagnostics[,3],cor((GY1all-GY1call),GXc0all))

  
  #We use GX and GGX constructed as instruments
  est1.2.1.2<-ivreg(Y1all ~ Xall + GY1call | Xall +  GXc0all + G2Xc0all)
  sest1.2.1.2<-summary(est1.2.1.2, diagnostic=TRUE)
  lest1.2.1.2 <- c(sest1.2.1.2$coefficients[,1],sest1.2.1.2$diagnostics[,3],cor((GY1all-GY1call),cbind(GXc0all,G2Xcall)))
  
  #different draw
  #We use GX constructed as instruments
  est1.2.2.1<-ivreg(Y1all ~ Xall + GY1call | Xall +  GXcall)
  sest1.2.2.1<-summary(est1.2.2.1, diagnostic=TRUE)
  lest1.2.2.1 <- c(sest1.2.2.1$coefficients[,1],sest1.2.2.1$diagnostics[,3],cor((GY1all-GY1call),GXcall))

  
  #We use GX and GGX constructed as instruments
  est1.2.2.2<-ivreg(Y1all ~ Xall + GY1call | Xall +  GXcall + G2Xcall)
  sest1.2.2.2<-summary(est1.2.2.2, diagnostic=TRUE)
  lest1.2.2.2 <- c(sest1.2.2.2$coefficients[,1],sest1.2.2.2$diagnostics[,3],cor((GY1all-GY1call),cbind(GXcall,G2Xcall)))
  
  # estimation section 3.2
  #GY is observed
  est2.1 <-
    ivreg(Y2all ~ Xall + GXall + GY2all |
            Xall  + GXall + G2Xcall)
  sest2.1 <- summary(est2.1, diagnostic = TRUE)
  lest2.1 <- c(sest2.1$coefficients[, 1], sest2.1$diagnostics[, 3])

  
  #GY is not observed
  est2.2 <-
    ivreg(Y2all ~ Xall + GXall +  GXc0all + GY2call |
            Xall  + GXall + GXc0all + G2Xcall)
  sest2.2 <- summary(est2.2, diagnostic = TRUE)
  lest2.2 <- c(sest2.2$coefficients[, 1], sest2.2$diagnostics[, 3])
  
  # estimation section 3.3
  est3.1 <-
    ivreg(Y3allm0 ~ Xallm0 + GXallm0 + GY3all |
            Xallm0 + GXallm0 + G2Xcallm0)
  sest3.1 <- summary(est3.1, diagnostic = TRUE)
  lest3.1 <- c(sest3.1$coefficients[, 1], sest3.1$diagnostics[, 3])
  
  #Y is not observed
  est3.2 <-
    ivreg(
      Y3allm0 ~ Xallm0 + GXallm0 + GXc0allm0 + GY3callm0 |
        Xallm0 + GXallm0 + GXc0allm0 + G2Xcallm0
    )
  sest3.2 <- summary(est3.2, diagnostic = TRUE)
  lest3.2 <- c(sest3.2$coefficients[, 1], sest3.2$diagnostics[, 3])
  
  sink(paste0("Aristide/Progress/Progress_", mpi.comm.rank(), "_", l, ".txt"))
  print(paste("Iteration:",l,"Model without contextual effects, Gy observed and GX is instrument"))
  print(lest1.1.1)
  print(paste("Iteration:",l,"Model without contextual effects, Gy observed GX, GGX are instruments"))
  print(lest1.1.2)
  print(paste("Iteration:",l,"Model without contextual effects, Gy not observed, same draws and GX is instrument"))
  print(lest1.2.1.1)
  print(paste("Iteration:",l,"Model without contextual effects, Gy not observed, same draws and GX, GGX are instruments"))
  print(lest1.2.1.2)
  print(paste("Iteration:",l,"Model without contextual effects, Gy not observed, different draws and GX is instrument"))
  print(lest1.2.2.1)
  print(paste("Iteration:",l,"Model without contextual effects, Gy not observed, different draws and GX, GGX are instruments"))
  print(lest1.2.2.2)
  print(paste("Iteration:", l, "Model with contextual effects, Gy is observed"))
  print(lest2.1)
  print(paste("Iteration:", l, "Model with contextual effects, Gy is non observed"))
  print(lest2.2)
  print(paste("Iteration:", l, "Model fixed effects, Gy is observed"))
  print(lest3.1)
  print(paste("Iteration:", l, "Model fixed effects, Gy is not observed"))
  print(lest3.2)
  sink()
  
  c(
    lest1.1.1,
    lest1.1.2,
    lest1.2.1.1,
    lest1.2.1.2,
    lest1.2.2.1,
    lest1.2.2.2,
    lest2.1,
    lest2.2,
    lest3.1,
    lest3.2
  )
}


fMC       <- function(iteration, parlambda) {
  
  ncores  <- mpi.universe.size() 
  print(ncores)
  mpi.spawn.Rslaves(nslaves = (ncores - 1)) 
  mpi.setup.rngstream() # generate random numbers
  
  mpi.bcast.cmd(library(AER))
  mpi.bcast.cmd(library(nuclearARD))
  mpi.bcast.cmd(library(doParallel))
  mpi.bcast.cmd(library(PartialNetwork))
  
  mpi.bcast.Robj2slave(obj,comm = 1,all = TRUE)
  
  if (dir.exists("Aristide/Progress")) {
    system("rm -r Aristide/Progress")
  }
  dir.create("Aristide/Progress")
  
  out.mc     <- mpi.parLapply(1:iteration, function(l)
    f(l, parlambda)) 
  
  saveRDS(out.mc, file = paste0("Save.nuclearARD.MC.l=", parlambda, ".RData"))
  mpi.close.Rslaves()
  mpi.quit()
}



fMC(1000, 200)
# fMC(1000, 600)
fMC(1000, 1374)





###########################################
rm(list = ls())

my.sum <- function(x) {
  out <- c(mean(x, na.rm = TRUE),
           sd(x, na.rm = TRUE),
           quantile(x, 0.25, na.rm = TRUE),
           median(x, na.rm = TRUE),
           quantile(x, 0.75, na.rm = TRUE))
  names(out) <- c("Mean", "Sd.", "1st Qu.", "Median", "3rd Qu.")
  return(out)
}


parlambda = 200
load(paste0("Save.nuclearARD.MC.l=", parlambda, ".RData"))

##Y is observe
### GX
(sum.lest1.1.1 <- t(as.matrix(apply(lest1.1.1[1:1000,] , 2, my.sum))))
### GX GGX
(sum.lest1.1.2 <- t(as.matrix(apply(lest1.1.2[1:1000,] , 2, my.sum))))
## y is not observed. Same draw
### GX
(sum.lest1.2.1.1 <- t(as.matrix(apply(lest1.2.1.1[1:1000,], 2, my.sum))))
### GX GGX
(sum.lest1.2.1.2 <- t(as.matrix(apply(lest1.2.1.2[1:1000,], 2, my.sum))))
## y is not observed. different draws
### GX
(sum.lest1.2.2.1 <- t(as.matrix(apply(lest1.2.2.1[1:1000,], 2, my.sum))))
### GX GGX
(sum.lest1.2.2.2 <- t(as.matrix(apply(lest1.2.2.2[1:1000,], 2, my.sum))))
## contextual effects
### y is observed
(sum.lest2.1 <- t(as.matrix(apply(lest2.1[1:1000,] , 2, my.sum))))
### y is not observed
(sum.lest2.2 <- t(as.matrix(apply(lest2.2[1:1000,] , 2, my.sum))))
## fixed effects
### y is observed
(sum.lest3.1 <- t(as.matrix(apply(lest3.1[1:1000,] , 2, my.sum))))
### y is not observed
(sum.lest3.2 <- t(as.matrix(apply(lest3.2[1:1000,] , 2, my.sum))))


write.csv(sum.lest1.1.1, file = "ard.nuc.sum.lest1.1.1.csv")
write.csv(sum.lest1.1.2, file = "ard.nuc.sum.lest1.1.2.csv")
#write.csv(sum.lest1.2.1.1, file = "ard.nuc.sum.lest1.2.1.1.csv")
#write.csv(sum.lest1.2.1.2, file = "ard.nuc.sum.lest1.2.1.2.csv")
write.csv(sum.lest1.2.2.1, file = "ard.nuc.sum.lest1.2.2.1.csv")
write.csv(sum.lest1.2.2.2, file = "ard.nuc.sum.lest1.2.2.2.csv")
write.csv(sum.lest2.1, file = "ard.nuc.sum.lest2.1.csv")
write.csv(sum.lest2.2, file = "ard.nuc.sum.lest2.2.csv")
write.csv(sum.lest3.1, file = "ard.nuc.sum.lest3.1.csv")
write.csv(sum.lest3.2, file = "ard.nuc.sum.lest3.2.csv")

