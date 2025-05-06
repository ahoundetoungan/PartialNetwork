library(PartialNetwork)
library(CDatanet)
library(doParallel)
library(dplyr)
library(ggplot2)
library(openxlsx)
rm(list = ls())

## Work directory 
## possible root 
proot <- c("~/Dropbox/Academy/1.Papers/Partial Network/Simulations/Monte Carlo",
           "~/PartialNetwork/Monte Carlo")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)


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
beta     <- c(3.8057, -0.0723, 0.1325)
gamma    <- c(0.0861, -0.0034)
alpha    <- 0.5378
rho      <- c(-2.348612, -0.699906, 0.403816)
se       <- 0.707

fMC      <- function(pobs){# pobs is the proportion of links observed
  # pobs   <- 0.25
  # Matrix X
  X        <- cbind(sample(agevalue, sum(N), replace = TRUE, prob = ageprob),
                    sample(femvalue, sum(N), replace = TRUE, prob = femprob))
  
  # Matrix X centered by group
  Xc       <- do.call(rbind, lapply(1:M, function(x){
    tmp    <- colMeans(X[c(CUMN[x] + 1):CUMN[x+1],])
    t(apply(X[c(CUMN[x] + 1):CUMN[x+1],], 1, function(w) w - tmp))}))
  
  # Simulate fixed effects as the median of X1
  eff      <- unlist(lapply(1:M, function(x) quantile(X[c(CUMN[x] + 1):CUMN[x+1],1], probs = 0.50)))
  eff      <- eff - mean(eff) + beta[1] # The average fixed effect = the intercept. This makes DGP with and without fixed effects comparable.
  eff      <- unlist(lapply(1:M, function(x) rep(eff[x], each = N[x])))
  
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
                     theta = c(alpha, 1, beta[-1], gamma, se))
  Gyf      <- yf$Gy
  yf       <- yf$y
  
  # Observed network
  nNet      <- nrow(Xnet) # network formation model sample size
  Aobs      <- sort(sample(1:nNet, round((1 - pobs)*nNet))) #Observed entries
  AOb       <- rep(0, nNet)
  AOb[Aobs] <- ynet[Aobs] 
  AOb       <- vec.to.mat(AOb, N, normalise = FALSE)
  GOb       <- norm.network(AOb)
  GOby      <- peer.avg(GOb, y)
  GObyf     <- peer.avg(GOb, yf)
  GObX      <- peer.avg(GOb, X)
  GObXc     <- peer.avg(GOb, Xc)
  GGObX     <- peer.avg(GOb, GObX)
  GGObXc    <- peer.avg(GOb, GObXc)
  
  # Put y, yf, Gy, Gyf, GX, GOby, GObyf, and GObX in the same dataset
  dataset           <- as.data.frame(cbind(y, yf, X, Gy, Gyf, GX, GOby, GObyf, GObX)) 
  colnames(dataset) <- c("y", "yf","X1","X2", "Gy", "Gyf", "GX1", "GX2", "GOby", "GObyf", "GObX1", "GObX2") 
  
  # Estimation the peer effect model using the observed network where Gy and GX are observed 
  # Without fixed effects
  W           <- solve(crossprod(cbind(1, X, GX, GGObX))/sum(N))
  gmm1        <- smmSAR(y ~ X1 + X2|Gy|GX1 + GX2, contextual = TRUE, dnetwork = AOb, cond.var = FALSE,
                        W = W, smm.ctr  = list(R = 1, print = FALSE), data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmm1["Gy"]  <- ifelse(gmm1["Gy"] < -1, -1, ifelse(gmm1["Gy"] > 1, 1, gmm1["Gy"]))
  names(gmm1) <- paste0("gmm-GyGXobs", names(gmm1))
  
  # With fixed effects
  Wf          <- solve(crossprod(cbind(Xc, GXc, GGObXc))/sum(N))
  gmmf1       <- smmSAR(yf ~ X1 + X2|Gyf|GX1 + GX2, contextual = TRUE, dnetwork = AOb, cond.var = FALSE,
                        fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 1, print = FALSE), 
                        data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmmf1["Gy"]  <- ifelse(gmmf1["Gy"] < -1, -1, ifelse(gmmf1["Gy"] > 1, 1, gmmf1["Gy"]))
  names(gmmf1) <- paste0("gmmf-GyGXobs", names(gmmf1))
  
  # Estimation the peer effect model using the observed network where Gy is observed and GX is not
  # Without fixed effects
  W           <- solve(crossprod(cbind(1, X, GObX, GGObX))/sum(N))
  gmm2        <- smmSAR(y ~ X1 + X2|Gy|GObX1 + GObX2, contextual = TRUE, dnetwork = AOb, cond.var = FALSE,
                        W = W, smm.ctr  = list(R = 1, print = FALSE), data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmm2["Gy"]  <- ifelse(gmm2["Gy"] < -1, -1, ifelse(gmm2["Gy"] > 1, 1, gmm2["Gy"]))
  names(gmm2) <- paste0("gmm-Gyobs.GXunobs", names(gmm2))
  
  # With fixed effects
  Wf          <- solve(crossprod(cbind(Xc, GObXc, GGObXc))/sum(N))
  gmmf2       <- smmSAR(yf ~ X1 + X2|Gyf|GObX1 + GObX2, contextual = TRUE, dnetwork = AOb, cond.var = FALSE,
                        fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 1, print = FALSE), 
                        data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmmf2["Gy"]  <- ifelse(gmmf2["Gy"] < -1, -1, ifelse(gmmf2["Gy"] > 1, 1, gmmf2["Gy"]))
  names(gmmf2) <- paste0("gmmf-Gyobs.GXunobs", names(gmmf2))
  
  
  # Estimation the peer effect model using the observed network where Gy is not observed and GX is
  # Without fixed effects
  W           <- solve(crossprod(cbind(1, X, GX, GGObX))/sum(N))
  gmm3        <- smmSAR(y ~ X1 + X2|GOby|GX1 + GX2, contextual = TRUE, dnetwork = AOb, cond.var = FALSE,
                        W = W, smm.ctr  = list(R = 1, print = FALSE), data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmm3["Gy"]  <- ifelse(gmm3["Gy"] < -1, -1, ifelse(gmm3["Gy"] > 1, 1, gmm3["Gy"]))
  names(gmm3) <- paste0("gmm-Gyunobs.GXobs", names(gmm3))
  
  # With fixed effects
  Wf          <- solve(crossprod(cbind(Xc, GXc, GGObXc))/sum(N))
  gmmf3       <- smmSAR(yf ~ X1 + X2|GObyf|GX1 + GX2, contextual = TRUE, dnetwork = AOb, cond.var = FALSE,
                        fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 1, print = FALSE), 
                        data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmmf3["Gy"]  <- ifelse(gmmf3["Gy"] < -1, -1, ifelse(gmmf3["Gy"] > 1, 1, gmmf3["Gy"]))
  names(gmmf3) <- paste0("gmmf-Gyunobs.GXobs", names(gmmf3))
  
  # Estimation the peer effect model using the observed network (Gy and GX are not observed)
  # Without fixed effects
  W           <- solve(crossprod(cbind(1, X, GObX, GGObX))/sum(N))
  gmm4        <- smmSAR(y ~ X1 + X2|GOby|GObX1 + GObX2, contextual = TRUE, dnetwork = AOb, cond.var = FALSE,
                        W = W, smm.ctr  = list(R = 1, print = FALSE), data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmm4["Gy"]  <- ifelse(gmm4["Gy"] < -1, -1, ifelse(gmm4["Gy"] > 1, 1, gmm4["Gy"]))
  names(gmm4) <- paste0("gmm-GyGXunobs", names(gmm4))
  
  # With fixed effects
  Wf          <- solve(crossprod(cbind(Xc, GObXc, GGObXc))/sum(N))
  gmmf4       <- smmSAR(yf ~ X1 + X2|GObyf|GObX1 + GObX2, contextual = TRUE, dnetwork = AOb, cond.var = FALSE,
                        fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 1, print = FALSE), 
                        data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmmf4["Gy"]  <- ifelse(gmmf4["Gy"] < -1, -1, ifelse(gmmf4["Gy"] > 1, 1, gmmf4["Gy"]))
  names(gmmf4) <- paste0("gmmf-GyGXunobs", names(gmmf4))
  
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
  smm1       <- smmSAR(y ~ X1 + X2|Gy|GX1 + GX2, contextual = TRUE, dnetwork = d.logit, cond.var = FALSE, 
                       W = W, smm.ctr  = list(R = 100, print = FALSE), data = dataset)$estimates
  smm1f      <- smmSAR(yf ~ X1 + X2|Gyf|GX1 + GX2, contextual = TRUE, dnetwork = d.logit, cond.var = FALSE, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = FALSE), 
                       data = dataset)$estimates
  # These estimators have closed-form expressions and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  smm1["Gy"]   <- ifelse(smm1["Gy"] < -1, -1, ifelse(smm1["Gy"] > 1, 1, smm1["Gy"]))
  smm1f["Gy"]  <- ifelse(smm1f["Gy"] < -1, -1, ifelse(smm1f["Gy"] > 1, 1, smm1f["Gy"]))
  names(smm1)  <- paste0("smm1", names(smm1))
  names(smm1f) <- paste0("smm1f", names(smm1f))
  
  # SMM with Gy observed and GX unobserved
  smm2       <- smmSAR(y ~ X1 + X2|Gy, contextual = TRUE, dnetwork = d.logit, cond.var = FALSE, 
                       W = W, smm.ctr  = list(R = 100, print = FALSE), data = dataset)$estimates
  smm2f      <- smmSAR(yf ~ X1 + X2|Gyf, contextual = TRUE, dnetwork = d.logit, cond.var = FALSE, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = FALSE), 
                       data = dataset)$estimates
  # These estimators have closed-form expressions and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  smm2["Gy"]   <- ifelse(smm2["Gy"] < -1, -1, ifelse(smm2["Gy"] > 1, 1, smm2["Gy"]))
  smm2f["Gy"]  <- ifelse(smm2f["Gy"] < -1, -1, ifelse(smm2f["Gy"] > 1, 1, smm2f["Gy"]))
  names(smm2)  <- paste0("smm2", names(smm2))
  names(smm2f) <- paste0("smm2f", names(smm2f))
  
  # SMM with Gy unobserved and GX observed
  smm3       <- smmSAR(y ~ X1 + X2||GX1 + GX2, contextual = TRUE, dnetwork = d.logit, cond.var = FALSE, 
                       W = W, smm.ctr  = list(R = 100, print = FALSE), data = dataset)$estimates
  smm3f      <- smmSAR(yf ~ X1 + X2||GX1 + GX2, contextual = TRUE, dnetwork = d.logit, cond.var = FALSE, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = FALSE), 
                       data = dataset)$estimates
  names(smm3)  <- paste0("smm3", names(smm3))
  names(smm3f) <- paste0("smm3f", names(smm3f))
  
  # SMM with Gy and GX unobserved
  smm4       <- smmSAR(y ~ X1 + X2, contextual = TRUE, dnetwork = d.logit, cond.var = FALSE, 
                       W = W, smm.ctr  = list(R = 100, print = FALSE), data = dataset)$estimates
  smm4f      <- smmSAR(yf ~ X1 + X2, contextual = TRUE, dnetwork = d.logit, cond.var = FALSE, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = FALSE), 
                       data = dataset)$estimates
  names(smm4)  <- paste0("smm4", names(smm4))
  names(smm4f) <- paste0("smm4f", names(smm4f))
  
  c(gmm1, gmmf1, gmm2, gmmf2, gmm3, gmmf3, gmm4, gmmf4, smm1, smm1f, smm2, smm2f, smm3, smm3f, smm4, smm4f)
}

fsimu  <- function(pobs, mc.cores){
  out  <- do.call(rbind, mclapply(1:1e3, function(i) fMC(pobs), mc.cores = mc.cores))
  saveRDS(out, file = paste0("Results/logit.mlink.pobs=", pobs, ".RDS"))
}
set.seed(123)
fsimu(0, 15)
fsimu(0.25, 15)
fsimu(0.50, 15)
fsimu(0.75, 15)

## Export result toward Excel
out <- cbind(t(apply(readRDS("Results/logit.mlink.pobs=0.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.mlink.pobs=0.25.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.mlink.pobs=0.5.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.mlink.pobs=0.75.RDS"), 2, function(x) c(mean(x), sd(x)))))
colnames(out) <- paste0(rep(c("mean.pobs.", "sd.pobs."), 4),
                        rep(c("0", "25", "50", "75"), each = 2))

## Table for the model without FE
outNOFE       <- out[c(sapply(c(12, seq(34, 78, 11)), function(x) x + 0:5)),
                     paste0(rep(c("mean.pobs.", "sd.pobs."), 2), rep(c(0, 25, 50, 75), each = 2))]
outNOFE       <- data.frame(Var = rownames(outNOFE), outNOFE)
outNOFE$Var   <- rep(paste0(c("$\\alpha = ", "$c = ", "$\\beta_1 = ", "$\\beta_2 = ", "$\\gamma_1 = ", "$\\gamma_2 = "), 
                            sprintf("%.3f", c(alpha, beta, gamma)), "$"), 6)
tpV           <- as_tibble(matrix(NA, 1, ncol(outNOFE)))
colnames(tpV) <- colnames(outNOFE)
tp            <- c("Classical IV: $\\mathbf{Gy}$ observed and $\\mathbf{GX}$ unobserved; Instruments: $\\mathbf{GX}^2$",
                   "Classical IV: $\\mathbf{Gy}$ and $\\mathbf{GX}$ unobserved; Instruments: $\\mathbf{GX}^2$",
                   "SGMM: $\\mathbf{Gy}$ and $\\mathbf{GX}$ observed; $T = 100$",
                   "SGMM: $\\mathbf{Gy}$ observed and $\\mathbf{GX}$ unobserved; $S = T = 100$",
                   "SGMM: $\\mathbf{Gy}$ unobserved and $\\mathbf{GX}$ observed; $S = T = 100$",
                   "SGMM: $\\mathbf{Gy}$ and $\\mathbf{GX}$ unobserved; $R = S = T = 100$")
for (s in 6:1) {
  tpV$Var     <- tp[s]
  outNOFE     <- outNOFE %>% add_row(tpV, .before = 1 + 6*(s - 1))
}
outNOFE       <- outNOFE %>% mutate(across(all_of(paste0("mean.pobs.", c(0, 25, 50, 75))),
                                           ~ ifelse(is.na(.), NA, sprintf("%.3f", .))),
                                    across(all_of(paste0("sd.pobs.", c(0, 25, 50, 75))),
                                           ~ ifelse(is.na(.), ., paste0("(", sprintf("%.3f", .), ")"))))

## Table for the model FE
outFE         <- out[c(sapply(c(12, seq(34, 78, 11)), function(x) x + 6:10)),
                     paste0(rep(c("mean.pobs.", "sd.pobs."), 2), rep(c(0, 25, 50, 75), each = 2))]
outFE         <- data.frame(Var = rownames(outFE), outFE)
outFE$Var     <- rep(paste0(c("$\\alpha = ", "$\\beta_1 = ", "$\\beta_2 = ", "$\\gamma_1 = ", "$\\gamma_2 = "), 
                            sprintf("%.3f", c(alpha, beta[-1], gamma)), "$"), 6)
tpV           <- as_tibble(matrix(NA, 1, ncol(outFE)))
colnames(tpV) <- colnames(outFE)
tp            <- c("Classical IV: $\\mathbf{Gy}$ observed and $\\mathbf{GX}$ unobserved; Instruments: $\\mathbf{GX}^2$",
                   "Classical IV: $\\mathbf{Gy}$ and $\\mathbf{GX}$ unobserved; Instruments: $\\mathbf{GX}^2$",
                   "SGMM: $\\mathbf{Gy}$ and $\\mathbf{GX}$ observed; $T = 100$",
                   "SGMM: $\\mathbf{Gy}$ observed and $\\mathbf{GX}$ unobserved; $S = T = 100$",
                   "SGMM: $\\mathbf{Gy}$ unobserved and $\\mathbf{GX}$ observed; $S = T = 100$",
                   "SGMM: $\\mathbf{Gy}$ and $\\mathbf{GX}$ unobserved; $R = S = T = 100$")
for (s in 6:1) {
  tpV$Var     <- tp[s]
  outFE       <- outFE %>% add_row(tpV, .before = 1 + 5*(s - 1))
}
outFE         <- outFE %>% mutate(across(all_of(paste0("mean.pobs.", c(0, 25, 50, 75))),
                                         ~ ifelse(is.na(.), NA, sprintf("%.3f", .))),
                                  across(all_of(paste0("sd.pobs.", c(0, 25, 50, 75))),
                                         ~ ifelse(is.na(.), NA, paste0("(", sprintf("%.3f", .), ")"))))
# Create a new workbook
wb            <- createWorkbook()
# Add data
addWorksheet(wb, "FullResult"); writeData(wb, sheet = "FullResult", x = data.frame(Var = rownames(out), out))
addWorksheet(wb, "NOFE"); writeData(wb, sheet = "NOFE", x = outNOFE)
addWorksheet(wb, "FE"); writeData(wb, sheet = "FE", x = outFE)
saveWorkbook(wb, "Results/logit.mlink.xlsx", overwrite = TRUE)

# Graph
mobs <- c(0, 0.25, 0.5, 0.75)
vmis <- c("0%", "25%", "50%", "75%")
data <- do.call(rbind, lapply(1:length(mobs), function(k) {
  est  <- apply(readRDS(paste0("Results/logit.mlink.pobs=", mobs[k], ".RDS")),
                2, function(x) c(mean(x), quantile(x, prob = c(0.025, 0.975))))
  data.frame(pmis  = vmis[k],
             spec  = rep(1:8, each = 2),
             model = c(rep("Classical IV: Gy, GX observed", 2),
                       rep("Classical IV: Gy observed, GX unobserved", 2),
                       rep("Classical IV: Gy, Gy unobserved, GX observed", 2),
                       rep("Classical IV: Gy, GX unobserved", 2),
                       rep("SGMM: Gy, GX observed", 2),
                       rep("SGMM: Gy observed, GX unobserved", 2),
                       rep("SGMM: Gy unobserved, GX observed", 2),
                       rep("SGMM: Gy, GX unobserved", 2)),
             FE    = rep(c(FALSE, TRUE), 8),
             coef  = est[1, c("gmm-GyGXobsGy", "gmmf-GyGXobsGy", 
                              "gmm-Gyobs.GXunobsGy", "gmmf-Gyobs.GXunobsGy", 
                              "gmm-Gyunobs.GXobsGy", "gmmf-Gyunobs.GXobsGy", 
                              "gmm-GyGXunobsGy", "gmmf-GyGXunobsGy", 
                              "smm1Gy", "smm1fGy", 
                              "smm2Gy", "smm2fGy", 
                              "smm3Gy", "smm3fGy", 
                              "smm4Gy", "smm4fGy")],
             IC1   = est[2, c("gmm-GyGXobsGy", "gmmf-GyGXobsGy", 
                              "gmm-Gyobs.GXunobsGy", "gmmf-Gyobs.GXunobsGy", 
                              "gmm-Gyunobs.GXobsGy", "gmmf-Gyunobs.GXobsGy", 
                              "gmm-GyGXunobsGy", "gmmf-GyGXunobsGy", 
                              "smm1Gy", "smm1fGy", 
                              "smm2Gy", "smm2fGy", 
                              "smm3Gy", "smm3fGy", 
                              "smm4Gy", "smm4fGy")],
             IC2   = est[3, c("gmm-GyGXobsGy", "gmmf-GyGXobsGy", 
                              "gmm-Gyobs.GXunobsGy", "gmmf-Gyobs.GXunobsGy", 
                              "gmm-Gyunobs.GXobsGy", "gmmf-Gyunobs.GXobsGy", 
                              "gmm-GyGXunobsGy", "gmmf-GyGXunobsGy", 
                              "smm1Gy", "smm1fGy", 
                              "smm2Gy", "smm2fGy", 
                              "smm3Gy", "smm3fGy", 
                              "smm4Gy", "smm4fGy")])
}))

data   <- data %>% mutate(Model = factor(spec, labels = unique(model)),
                          Which = "Missing links")


(NOFE  <- ggplot(data %>% filter(FE == FALSE, spec %in% c(4, 8)), aes(x = pmis, colour = Model)) +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "gray") +
    geom_errorbar(width=.2, aes(ymin = IC1, ymax = IC2),
                  position = position_dodge(width = 0.3)) +
    geom_point(aes(y = coef, shape = Model), position = position_dodge(width = 0.3)) +
    theme_bw() +
    facet_wrap(~Which) +
    annotate("text", x = 0.42, y = alpha, label = bquote(alpha[0] == .(round(alpha, 3))), 
             hjust = 0, vjust = -0.2, size = 2.8, color = "#333") +
    xlab("Proportion of missing links") + ylab("Peer effect estimate") +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)))
ggsave("mc_mlinkNOFE.pdf", plot = NOFE, device = "pdf", width = 6, height = 3.5)

# Save Data to merge with Missclassified links data
saveRDS(data, file = "Results/Missinggraph.RDS")

(WIFE  <- ggplot(data %>% filter(FE == TRUE, spec %in% c(4, 8)), aes(x = pmis, colour = Model)) +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "gray") +
    geom_errorbar(width=.2, aes(ymin = IC1, ymax = IC2),
                  position = position_dodge(width = 0.3)) +
    geom_point(aes(y = coef, shape = Model), position = position_dodge(width = 0.3)) +
    theme_bw() +
    facet_wrap(~Which) +
    annotate("text", x = 0.42, y = alpha, label = bquote(alpha[0] == .(round(alpha, 3))), 
             hjust = 0, vjust = -0.2, size = 2.8, color = "#333") +
    xlab("Proportion of missing links") + ylab("Peer effect estimate") +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)))
