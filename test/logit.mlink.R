library(PartialNetwork)
library(CDatanet)
library(doParallel)

setwd("~/Dropbox/Academy/1.Papers/Partial Network/Simulations/Monte Carlo")

rm(list = ls())

M        <- 100           # Number of groups
N        <- rep(30,M)   # Group size
CUMN     <- c(0, cumsum(N))

# Distribution of age
agevalue <- 10:19
ageprob  <- c(0.0006397953, 0.0102367242, 0.2271273193, 0.3458093410, 0.1896992962, 0.0854126679,
              0.0713371721, 0.0550223928, 0.0131158029, 0.0015994882) #distribution

# Distribution of female
femvalue <- 0:1
femprob  <- c(0.4596929, 0.5403071)

# Parameters
beta     <- c(1.2733, -0.068, 0.122)
gamma    <- c(0.103, -0.010)
alpha    <- 0.683
rho      <- c(-2.346986, 0.404115, -0.699205)
se       <- 0.7066

fMC      <- function(pobs){# pobs is the proportion of links observed
  # pobs   <- 0.25
  # Matrix X
  X        <- cbind(sample(agevalue, sum(N), replace = TRUE, prob = ageprob),
                    sample(femvalue, sum(N), replace = TRUE, prob = femprob))
  
  # Matrix X centered by group
  Xc       <- do.call(rbind, lapply(1:M, function(x){
    tmp    <- colMeans(X[c(CUMN[x] + 1):CUMN[x+1],])
    t(apply(X[c(CUMN[x] + 1):CUMN[x+1],], 1, function(w) w - tmp))}))
  
  # Simulate fixed effect as 25th centile of X2
  eff      <- unlist(lapply(1:M, function(x)
    rep(quantile(X[c(CUMN[x] + 1):CUMN[x+1],1], probs = 0.25), each = N[x])))
  
  # True network distribution
  X1l      <- lapply(1:M, function(x) X[c(CUMN[x] + 1):CUMN[x+1],1])
  X2l      <- lapply(1:M, function(x) X[c(CUMN[x] + 1):CUMN[x+1],2])
  X1.mat   <- lapply(1:M, function(m) {
    matrix(kronecker(X1l[[m]], X1l[[m]], FUN = function(x, y) abs(x - y)), N[m])})
  X2.mat   <- lapply(1:M, function(m) {
    matrix(kronecker(X2l[[m]], X2l[[m]], FUN = function(x, y) as.integer(x == y)), N[m])})
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
  y        <- simsar(~ X + GX, Glist = G0norm, theta = c(alpha, beta, gamma, se))
  Gy       <- y$Gy
  y        <- y$y
  
  # Simulate y with fixed effects
  yf       <- simsar(~-1 + eff + X + GX, Glist = G0norm, 
                     theta = c(alpha, 0.5, beta[-1], gamma, se))
  Gyf      <- yf$Gy
  yf       <- yf$y
  
  # Put y, yf, Gy, Gyf, GX in the same dataset
  dataset           <- as.data.frame(cbind(y, yf, X, Gy, Gyf, GX)) 
  colnames(dataset) <- c("y", "yf","X1","X2", "Gy", "Gyf", "GX1", "GX2") 
  
  # Observed network
  nNet      <- nrow(Xnet) # network formation model sample size
  Aobs      <- sort(sample(1:nNet, round((1 - pobs)*nNet))) #Observed entries
  AOb       <- rep(0, nNet)
  AOb[Aobs] <- ynet[Aobs] 
  AOb       <- vec.to.mat(AOb, N, normalise = FALSE)
  GOb       <- norm.network(AOb)
  GObX      <- peer.avg(GOb, X)
  GObXc     <- peer.avg(GOb, Xc)
  GGObX     <- peer.avg(GOb, GObX)
  GGObXc    <- peer.avg(GOb, GObXc)
  
  # Estimation the peer effect model using the observed network 
  # Without fixed effects
  W           <- solve(crossprod(cbind(1, X, GObX, GGObX))/sum(N))
  gmm         <- smmSAR(y ~ X1 + X2, contextual = T, dnetwork = AOb, cond.var = FALSE,
                        W = W, smm.ctr  = list(R = 1, print = F), data = dataset)$estimates
  names(gmm)  <- paste0("gmm", names(gmm))
  
  # With fixed effects
  Wf          <- solve(crossprod(cbind(Xc, GObXc, GGObXc))/sum(N))
  gmmf        <- smmSAR(yf ~ X1 + X2, contextual = T, dnetwork = AOb, cond.var = FALSE,
                        fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 1, print = F), 
                        data = dataset)$estimates
  names(gmmf) <- paste0("gmmf", names(gmmf))
  
  
  # Estimation of the network distribution
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
  saveRDS(out, file = paste0("logit.mlink.pobs=", pobs, ".RDS"))
}
fsimu(0, 15)
fsimu(0.05, 15)
fsimu(0.10, 15)
fsimu(0.25, 15)
fsimu(0.50, 15)
fsimu(0.75, 15)

out <- cbind(t(apply(readRDS("logit.mlink.pobs=0.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("logit.mlink.pobs=0.05.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("logit.mlink.pobs=0.1.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("logit.mlink.pobs=0.25.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("logit.mlink.pobs=0.5.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("logit.mlink.pobs=0.75.RDS"), 2, function(x) c(mean(x), sd(x)))))

colnames(out) <- paste0(rep(c("mean.pobs=", "sd.pobs="), 5),
                        rep(c("0%", "5%", "10%", "25%", "50%", "75%"), each = 2))

write.csv(out, file = "Results/logit.mlink.csv")

# # Graph
# library(ggplot2)
# library(dplyr)
# est  <- apply(readRDS("../../../Simulations/Monte Carlo/logit.pobs=0.1.RDS"), 2, 
#               function(x) c(mean(x), quantile(x, prob = c(0.025, 0.975))))
# data <- data.frame(pmix  = rep("0%", 2),
#                    spec  = rep(0, 2),
#                    model = rep("Classical IV", 2),
#                    FE    = c(FALSE, TRUE),
#                    coef  = est[1, c("gmmGy", "gmmfGy")],
#                    IC1   = est[2, c("gmmGy", "gmmfGy")],
#                    IC2   = est[3, c("gmmGy", "gmmfGy")])
# mobs <- c(0.75, 0.5, 0.25, 0.1, 0.05)
# vmis <- c("25%", "50%", "75%", "90%", "95%")
# for (k in 1:length(mobs)) {
#   est  <- apply(readRDS(paste0("../../../Simulations/Monte Carlo/logit.pobs=", mobs[k], ".RDS")), 
#                 2, function(x) c(mean(x), quantile(x, prob = c(0.025, 0.975))))
#   data  <- data %>% 
#     bind_rows(data.frame(pmix  = rep(vmis[k], 4),
#                          spec  = rep(c(1, 3, 2, 4), each = 2),
#                          model = c(rep("SGMM: Gy, GX observed", 2),
#                                    rep("SGMM: Gy observed, GX unobserved", 2),
#                                    rep("SGMM: Gy unobserved, GX observed", 2),
#                                    rep("SGMM: Gy, GX unobserved", 2)),
#                          FE    = rep(c(FALSE, TRUE),  4),
#                          coef  = est[1, c("smm1Gy", "smm1fGy", "smm2Gy", "smm2fGy",
#                                           "smm3Gy", "smm3fGy", "smm4Gy", "smm4fGy")],
#                          IC1   = est[2, c("smm1Gy", "smm1fGy", "smm2Gy", "smm2fGy",
#                                           "smm3Gy", "smm3fGy", "smm4Gy", "smm4fGy")],
#                          IC2   = est[3, c("smm1Gy", "smm1fGy", "smm2Gy", "smm2fGy",
#                                           "smm3Gy", "smm3fGy", "smm4Gy", "smm4fGy")]))
# }
# 
# data   <- data %>% mutate(Model = factor(spec, labels = unique(model)))
# 
# ggplot(data %>% filter(pmix %in% c("0%", "25%", "50%", "75%"), 
#                        FE == FALSE, spec %in% c(0:4)), aes(x = pmix, colour = Model)) + 
#   geom_errorbar(width=.2, aes(ymin = IC1, ymax = IC2),
#                 position = position_dodge(width = 0.3)) +
#   geom_point(aes(y = coef, shape = Model), position = position_dodge(width = 0.3)) + 
#   theme_bw() +
#   xlab("Proportion of missing links") + ylab("Peer effect estimate") + 
#   theme(legend.title = element_blank(), legend.position = "bottom") + 
#   guides(colour = guide_legend(nrow = 2, byrow = FALSE))
# 
# ggplot(data %>% filter(pmix %in% c("0%", "25%", "50%", "75%"), 
#                        FE == TRUE, spec %in% c(0:4)), aes(x = pmix, colour = Model)) + 
#   geom_errorbar(width=.2, aes(ymin = IC1, ymax = IC2),
#                 position = position_dodge(width = 0.3)) +
#   geom_point(aes(y = coef, shape = Model), position = position_dodge(width = 0.3)) + 
#   theme_bw() +
#   xlab("Proportion of missing links") + ylab("Peer effect estimate") + 
#   theme(legend.title = element_blank(), legend.position = "bottom") + 
#   guides(colour = guide_legend(nrow = 2, byrow = FALSE))
# 


