library(PartialNetwork)
library(CDatanet)
library(doParallel)

rm(list = ls())

M        <- 100           # Number of groups
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

out <- cbind(t(apply(`logit.pobs=0.05`, 2, function(x) c(mean(x), sd(x)))),
             t(apply(`logit.pobs=0.1`, 2, function(x) c(mean(x), sd(x)))),
             t(apply(`logit.pobs=0.25`, 2, function(x) c(mean(x), sd(x)))),
             t(apply(`logit.pobs=0.5`, 2, function(x) c(mean(x), sd(x)))),
             t(apply(`logit.pobs=0.75`, 2, function(x) c(mean(x), sd(x)))))

colnames(out) <- paste0(rep(c("mean.pobs=", "sd.pobs="), 5),
                        rep(c("5%", "10%", "25%", "50%", "75%"), each = 2))

write.csv(out, file = "../../../Simulations/Monte Carlo/Results/logit.csv")

# Graph
library(ggplot2)
library(dplyr)
est  <- apply(readRDS("../../../Simulations/Monte Carlo/logit.pobs=0.1.RDS"), 2, 
              function(x) c(mean(x), quantile(x, prob = c(0.025, 0.975))))
data <- data.frame(pmix  = rep("0%", 2),
                   spec  = rep(0, 2),
                   model = rep("Classical IV", 2),
                   FE    = c(FALSE, TRUE),
                   coef  = est[1, c("gmmGy", "gmmfGy")],
                   IC1   = est[2, c("gmmGy", "gmmfGy")],
                   IC2   = est[3, c("gmmGy", "gmmfGy")])
mobs <- c(0.75, 0.5, 0.25, 0.1, 0.05)
vmis <- c("25%", "50%", "75%", "90%", "95%")
for (k in 1:length(mobs)) {
  est  <- apply(readRDS(paste0("../../../Simulations/Monte Carlo/logit.pobs=", mobs[k], ".RDS")), 
                2, function(x) c(mean(x), quantile(x, prob = c(0.025, 0.975))))
  data  <- data %>% 
    bind_rows(data.frame(pmix  = rep(vmis[k], 4),
                         spec  = rep(c(1, 3, 2, 4), each = 2),
                         model = c(rep("SGMM: Gy, GX observed", 2),
                                   rep("SGMM: Gy observed, GX unobserved", 2),
                                   rep("SGMM: Gy unobserved, GX observed", 2),
                                   rep("SGMM: Gy, GX unobserved", 2)),
                         FE    = rep(c(FALSE, TRUE),  4),
                         coef  = est[1, c("smm1Gy", "smm1fGy", "smm2Gy", "smm2fGy",
                                          "smm3Gy", "smm3fGy", "smm4Gy", "smm4fGy")],
                         IC1   = est[2, c("smm1Gy", "smm1fGy", "smm2Gy", "smm2fGy",
                                          "smm3Gy", "smm3fGy", "smm4Gy", "smm4fGy")],
                         IC2   = est[3, c("smm1Gy", "smm1fGy", "smm2Gy", "smm2fGy",
                                          "smm3Gy", "smm3fGy", "smm4Gy", "smm4fGy")]))
}

data   <- data %>% mutate(Model = factor(spec, labels = unique(model)))

ggplot(data %>% filter(pmix %in% c("0%", "25%", "50%", "75%"), 
                       FE == FALSE, spec %in% c(0:4)), aes(x = pmix, colour = Model)) + 
  geom_errorbar(width=.2, aes(ymin = IC1, ymax = IC2),
                position = position_dodge(width = 0.3)) +
  geom_point(aes(y = coef, shape = Model), position = position_dodge(width = 0.3)) + 
  theme_bw() +
  xlab("Proportion of missing links") + ylab("Peer effect estimate") + 
  theme(legend.title = element_blank(), legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 2, byrow = FALSE))

ggplot(data %>% filter(pmix %in% c("0%", "25%", "50%", "75%"), 
                       FE == TRUE, spec %in% c(0:4)), aes(x = pmix, colour = Model)) + 
  geom_errorbar(width=.2, aes(ymin = IC1, ymax = IC2),
                position = position_dodge(width = 0.3)) +
  geom_point(aes(y = coef, shape = Model), position = position_dodge(width = 0.3)) + 
  theme_bw() +
  xlab("Proportion of missing links") + ylab("Peer effect estimate") + 
  theme(legend.title = element_blank(), legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 2, byrow = FALSE))



