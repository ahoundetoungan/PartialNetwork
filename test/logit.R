library(PartialNetwork)
library(CDatanet)
library(doParallel)

rm(list = ls())

M        <- 200           # Number of groups
N        <- rep(30,M)   # Group size
CUMN     <- c(0, cumsum(N))

# Parameters
beta     <- c(2, 1, 1.5)
gamma    <- c(5, -3)
alpha    <- 0.4
rho      <- c(0.8, -0.2, 0.1)
se       <- 1

fMC      <- function(pobs){# pobs is the proportion of links observed
  # Matrix X
  X        <- cbind(rnorm(sum(N), 0, 5),rpois(sum(N), 6))
  
  # Matrix X centered by group
  Xc       <- do.call(rbind, lapply(1:M, function(x){
    tmp    <- colMeans(X[c(CUMN[x] + 1):CUMN[x+1],])
    t(apply(X[c(CUMN[x] + 1):CUMN[x+1],], 1, function(w) w - tmp))}))
  
  # Simulate fixed effect as 25th centile of X2
  eff      <- unlist(lapply(1:M, function(x)
    rep(quantile(X[c(CUMN[x] + 1):CUMN[x+1],2], probs = 0.25), each = N[x])))
  
  # True network distribution
  X1l      <- lapply(1:M, function(x) X[c(CUMN[x] + 1):CUMN[x+1],1])
  X2l      <- lapply(1:M, function(x) X[c(CUMN[x] + 1):CUMN[x+1],2])
  dist.net <- function(x, y) abs(x - y)
  X1.mat   <- lapply(1:M, function(m) {
    matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), N[m])})
  X2.mat   <- lapply(1:M, function(m) {
    matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), N[m])})
  Xnet     <- as.matrix(cbind("Const" = 1, 
                              "dX1"   = mat.to.vec(X1.mat),
                              "dX2"   = mat.to.vec(X2.mat)))
  ynet     <- Xnet %*% rho
  ynet     <- c(1*((ynet + rlogis(length(ynet))) > 0))
  G0       <- vec.to.mat(ynet, N, normalise = FALSE)
  G0norm   <- norm.network(G0)
  
  # GX, GGX also with the centered X
  GX       <- peer.avg(G0norm, X)
  GXc      <- peer.avg(G0norm, Xc)
  GGX      <- peer.avg(G0norm, GX)
  GGXc     <- peer.avg(G0norm, GXc)
  
  # Simulate y without fixed effects
  y        <- simsar(~ X, contextual = TRUE, Glist = G0norm, 
                     theta = c(alpha, beta, gamma, se))
  Gy       <- y$Gy
  y        <- y$y
  
  # Simulate y with fixed effects
  yf       <- simsar(~-1 + eff + X | X, Glist = G0norm, 
                     theta = c(alpha, 0.5, beta[-1], gamma, se))
  Gyf      <- yf$Gy
  yf       <- yf$y
  
  # Put y, yf, Gy, Gyf, GX in the same dataset
  dataset           <- as.data.frame(cbind(y, yf, X, Gy, Gyf, GX)) 
  colnames(dataset) <- c("y", "yf","X1","X2", "Gy", "Gyf", "GX1", "GX2") 
  
  # Estimation the peer effect model using the true network 
  # Without fixed effects
  W           <- solve(crossprod(cbind(1, X, GX, GGX))/sum(N))
  gmm         <- smmSAR(y ~ X1 + X2, contextual = T, dnetwork = G0, cond.var = FALSE,
                        W = W, smm.ctr  = list(R = 1, print = F), data = dataset)$estimates
  names(gmm)  <- paste0("gmm", names(gmm))
  
  # With fixed effects
  Wf          <- solve(crossprod(cbind(Xc, GXc, GGXc))/sum(N))
  gmmf        <- smmSAR(yf ~ X1 + X2, contextual = T, dnetwork = G0, cond.var = FALSE,
                        fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 1, print = F), 
                        data = dataset)$estimates
  names(gmmf) <- paste0("gmmf", names(gmmf))
  
  # Estimation of the network distribution
  nNet      <- nrow(Xnet) # network formation model sample size
  Aobs      <- sort(sample(1:nNet, round(pobs*nNet))) #Observed entries
  logestim  <- glm(ynet[Aobs] ~ -1 + Xnet[Aobs,], family = binomial(link = "logit"))
  slogestim <- summary(logestim)
  rho.est   <- logestim$coefficients
  hpl       <- c(1/(1 + exp(-as.matrix(Xnet)%*%rho.est)))
  hpl[Aobs] <- ynet[Aobs]
  d.logit   <- vec.to.mat(hpl, N)
  
  # Cumpute GX and GGX using simulated network 
  Gsim        <- sim.network(d.logit, normalise = TRUE)
  GsimX       <- peer.avg(Gsim, X)
  GsimGsimX   <- peer.avg(Gsim, GsimX)
  GsimXc      <- peer.avg(Gsim, Xc)
  GsimGsimXc  <- peer.avg(Gsim, GsimXc)
  W           <- solve(crossprod(cbind(1, X, GsimX, GsimGsimX))/sum(N))
  Wf          <- solve(crossprod(cbind(Xc, GsimXc, GsimGsimXc))/sum(N))
  
  # SMM with Gy and GX observed
  smm1       <- smmSAR(y ~ X1 + X2|Gy|GX1 + GX2, contextual = T, dnetwork = d.logit, 
                       W = W, smm.ctr  = list(R = 100, print = F), data = dataset)$estimates
  smm1f      <- smmSAR(yf ~ X1 + X2|Gyf|GX1 + GX2, contextual = T, dnetwork = d.logit, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = F), 
                       data = dataset)$estimates
  names(smm1)  <- paste0("smm1", names(smm1))
  names(smm1f) <- paste0("smm1f", names(smm1f))
  
  # SMM with Gy observed and GX unobserved
  smm2       <- smmSAR(y ~ X1 + X2|Gy, contextual = T, dnetwork = d.logit, 
                       W = W, smm.ctr  = list(R = 100, print = F), data = dataset)$estimates
  smm2f      <- smmSAR(yf ~ X1 + X2|Gyf, contextual = T, dnetwork = d.logit, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = F), 
                       data = dataset)$estimates
  names(smm2)  <- paste0("smm2", names(smm2))
  names(smm2f) <- paste0("smm2f", names(smm2f))
  
  # SMM with Gy unobserved and GX observed
  smm3       <- smmSAR(y ~ X1 + X2||GX1 + GX2, contextual = T, dnetwork = d.logit, 
                       W = W, smm.ctr  = list(R = 100, print = F), data = dataset)$estimates
  smm3f      <- smmSAR(yf ~ X1 + X2||GX1 + GX2, contextual = T, dnetwork = d.logit, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = F), 
                       data = dataset)$estimates
  names(smm3)  <- paste0("smm3", names(smm3))
  names(smm3f) <- paste0("smm3f", names(smm3f))
  
  # SMM with Gy and GX unobserved
  smm4       <- smmSAR(y ~ X1 + X2, contextual = T, dnetwork = d.logit, 
                       W = W, smm.ctr  = list(R = 100, print = F), data = dataset)$estimates
  smm4f      <- smmSAR(yf ~ X1 + X2, contextual = T, dnetwork = d.logit, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = F), 
                       data = dataset)$estimates
  names(smm4)  <- paste0("smm4", names(smm4))
  names(smm4f) <- paste0("smm4f", names(smm4f))
  
  c(gmm, gmmf, smm1, smm1f, smm2, smm2f, smm3, smm3f, smm4, smm4f)}

fsimu  <- function(pobs, mc.cores){
  out  <- do.call(rbind, mclapply(1:1e3, function(i) fMC(pobs), mc.cores = mc.cores))
  saveRDS(out, file = paste0("logit.pobs=", pobs, ".RDS"))
}

fsimu(0.05, 8)
fsimu(0.10, 8)
fsimu(0.25, 8)
fsimu(0.50, 8)
fsimu(0.75, 8)
