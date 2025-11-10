## ----setup, include=FALSE-----------------------------------------------------
rmarkdown::find_pandoc(version = '2.9.2.1')
knitr::opts_chunk$set(echo = TRUE)

## ----packages, include=FALSE, eval=TRUE---------------------------------------
# Check availability of packages 
has_ivreg <- requireNamespace("ivreg", quietly = TRUE) 
has_ggplt <- requireNamespace("ggplot2", quietly = TRUE) 

## ----IV1, echo = TRUE, eval = TRUE--------------------------------------------
library(PartialNetwork)
set.seed(123)
# Number of groups
M             <- 30
# size of each group
N             <- rep(50,M)
# individual effects
beta          <- c(2,1,1.5)
# endogenous effects
alpha         <- 0.4
# std-dev errors
se            <- 1
# network distribution
distr         <- runif(sum(N*(N-1)))
distr         <- vec.to.mat(distr, N, normalise = FALSE)
# covariates
X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# true network
G0            <- sim.network(distr)
# normalise
G0norm        <- norm.network(G0)
# simulate dependent variable
y             <- simSAR(~ X, Glist = G0norm, parms = c(alpha, beta), 
                        epsilon = rnorm(sum(N), sd = se))
y             <- y$y
# generate instruments
instr         <- sim.IV(distr, X, y, replication = 1, power = 1)
GY1c1         <- instr[[1]]$G1y       # proxy for Gy (draw 1)
GXc1          <- instr[[1]]$G1X[,,1]  # proxy for GX (draw 1)
GXc2          <- instr[[1]]$G2X[,,1]  # proxy for GX (draw 2)
# build dataset
# keep only instrument constructed using a different draw than the one used to proxy Gy
dataset           <- as.data.frame(cbind(y, X, GY1c1, GXc1, GXc2))
colnames(dataset) <- c("y","X1","X2","G1y", "G1X1", "G1X2", "G2X1", "G2X2")

## ----IV2, echo = TRUE, eval = has_ivreg, message=FALSE------------------------
library(ivreg)
# Same draws
out.iv1           <- ivreg(y ~ X1 + X2 + G1y | X1 + X2 + G1X1 + G1X2, data = dataset)
summary(out.iv1)

## ----IV3, echo = TRUE, eval = has_ivreg, message=FALSE------------------------
# Different draws
out.iv2           <- ivreg(y ~ X1 + X2 + G1y | X1 + X2 + G2X1 + G2X2, data = dataset)
summary(out.iv2)

## ----IV3bias------------------------------------------------------------------
ddS     <- as.matrix(cbind(1, dataset[,c("X1", "X2", "G1y")]))         #\ddot{S}
dZ      <- as.matrix(cbind(1, dataset[,c("X1", "X2", "G2X1", "G2X2")]))#\dot{Z}
dZddS   <- crossprod(dZ, ddS)/sum(N)                                  
W       <- solve(crossprod(dZ)/sum(N))                                 
matM    <- solve(crossprod(dZddS, W%*%dZddS), crossprod(dZddS, W))    
maxbias <- apply(sapply(1:1000, function(...){
  dddGy <- peer.avg(sim.network(distr, normalise = TRUE) , y)
  abs(matM%*%crossprod(dZ, dddGy - dataset$G1y)/sum(N))
}), 1, max); names(maxbias) <- c("(Intercept)", "X1", "X2", "G1y")
{cat("Maximal absolute bias\n"); print(maxbias)}

## ----IV4, echo = TRUE, eval = TRUE--------------------------------------------
set.seed(123)
# Number of groups
M             <- 30
# size of each group
N             <- rep(50,M)
# individual effects
beta          <- c(2,1,1.5)
# contextual effects
gamma         <- c(5, -3)
# endogenous effects
alpha         <- 0.4
# std-dev errors
se            <- 1
# network distribution
distr         <- runif(sum(N*(N-1)))
distr         <- vec.to.mat(distr, N, normalise = FALSE)
# covariates
X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# true network
G0            <- sim.network(distr)
# normalise
G0norm        <- norm.network(G0)
# GX
GX            <- peer.avg(G0norm, X)
# simulate dependent variable
y             <- simSAR(~ X + GX, Glist = G0norm, parms = c(alpha, beta, gamma), 
                        epsilon = rnorm(sum(N), sd = se))
y             <- y$y
# generate instruments
# we need power = 2 for models with contextual effetcs
instr         <- sim.IV(distr, X, y, replication = 1, power = 2)
GY1c1         <- instr[[1]]$G1y       # proxy for Gy (draw 1)
GXc1          <- instr[[1]]$G1X[,,1]  # proxy for GX (draw 1)
GXc2          <- instr[[1]]$G2X[,,1]  # proxy for GX (draw 2)
GXc2sq        <- instr[[1]]$G2X[,,2]  # proxy for G^2X (draw 2)
# build dataset
# keep only instrument constructed using a different draw than the one used to proxy Gy
dataset           <- as.data.frame(cbind(y, X, GX, GY1c1, GXc1, GXc2, GXc2sq))
colnames(dataset) <- c("y","X1","X2", "GX1", "GX2", "G1y", "G1X1", "G1X2", "G2X1", "G2X2",
                       "G2X1sq", "G2X2sq")

## ----IV5, echo = TRUE, eval = has_ivreg, message=FALSE------------------------
# Different draws
out.iv2           <- ivreg(y ~ X1 + X2 + GX1 + GX2 + G1X1 + G1X2 + G1y | X1 + X2 + GX1 +
                             GX2 + G1X1 + G1X2 + G2X1 + G2X2 + G2X1sq + G2X2sq,
                           data = dataset)
summary(out.iv2)

## ----IV5bias------------------------------------------------------------------
ddS     <- as.matrix(cbind(1, dataset[,c("X1", "X2", "GX1", "GX2", "G1X1", "G1X2",
                                         "G1y")]))                
dZ      <- as.matrix(cbind(1, dataset[,c("X1", "X2", "GX1", "GX2", "G1X1",
                                         "G1X2", "G2X1", "G2X2", "G2X1sq", "G2X2sq")]))
dZddS   <- crossprod(dZ, ddS)/sum(N)                                  
W       <- solve(crossprod(dZ)/sum(N))                                 
matM    <- solve(crossprod(dZddS, W%*%dZddS), crossprod(dZddS, W))    
maxbias <- apply(sapply(1:1000, function(...){
  dddGy <- peer.avg(sim.network(distr, normalise = TRUE) , y)
  abs(matM%*%crossprod(dZ, dddGy - dataset$G1y)/sum(N))
}), 1, max); names(maxbias) <- c("(Intercept)", "X1", "X2", "GX1", "GX2", "G1X1",
                                 "G1X2", "G1y")
{cat("Maximal absolute bias\n"); print(maxbias)}

## ----smm1, echo = TRUE, eval = TRUE-------------------------------------------
set.seed(123)
# Number of groups
M             <- 100
# size of each group
N             <- rep(30,M)
# individual effects
beta          <- c(2, 1, 1.5, 5, -3)
# endogenous effects
alpha         <- 0.4
# std-dev errors
se            <- 1
# network distribution
distr         <- runif(sum(N*(N-1)))
distr         <- vec.to.mat(distr, N, normalise = FALSE)
# covariates
X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# true network
G0            <- sim.network(distr)
# normalise
G0norm        <- norm.network(G0)
# Matrix GX
GX            <- peer.avg(G0norm, X)
# simulate dependent variable
y             <- simSAR(~ X + GX, Glist = G0norm, parms = c(alpha, beta), 
                        epsilon = rnorm(sum(N), sd = se))
Gy            <- y$Gy
y             <- y$y
# build dataset
dataset           <- as.data.frame(cbind(y, X, Gy, GX))
colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2")

## ----Smmload, echo = FALSE, eval = TRUE, message=FALSE------------------------
load(url('https://raw.githubusercontent.com/ahoundetoungan/PartialNetwork/master/datavignettes/smm.rda', "rb"))

## ----Smm2, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smm1       <- smmSAR(y ~ X1 + X2 | Gy | GX1 + GX2, dnetwork = distr, contextual = T,
#                          smm.ctr  = list(R = 1, print = F), data = dataset)
# summary(out.smm1)

## ----Smm2a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smm1)

## ----Smm3, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smm2       <- smmSAR(y ~ X1 + X2 || GX1 + GX2, dnetwork = distr, contextual = T,
#                          smm.ctr  = list(R = 1, print = F), data = dataset)
# summary(out.smm2)

## ----Smm3a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smm2)

## ----Smm4, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smm3       <- smmSAR(y ~ X1 + X2 | Gy, dnetwork = distr, contextual = T,
#                          smm.ctr  = list(R = 100, print = F), data = dataset)
# summary(out.smm3)

## ----Smm4a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smm3)

## ----Smm5, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smm4       <- smmSAR(y ~ X1 + X2, dnetwork = distr, contextual = T,
#                          smm.ctr  = list(R = 100, print = F), data = dataset)
# summary(out.smm4)

## ----Smm5a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smm4)

## ----smm6b, echo = TRUE, eval = TRUE------------------------------------------
set.seed(123)
# Number of groups
M             <- 200
# size of each group
N             <- rep(30,M)
# individual effects
beta          <- c(1, 1, 1.5, 5, -3)
# endogenous effects
alpha         <- 0.4
# std-dev errors
se            <- 1
# network distribution
distr         <- runif(sum(N*(N-1)))
distr         <- vec.to.mat(distr, N, normalise = FALSE)
# covariates
X             <- cbind(rnorm(sum(N),0,5), rpois(sum(N),7))
# Groups' fixed effects
# In order to have groups' heterogeneity correlated to X (fixed effects),
# We consider the quantile of X2 at 25% in the group
eff           <- unlist(lapply(1:M, function(x)
  rep(quantile(X[(c(0, cumsum(N))[x]+1):(cumsum(N)[x]),2], probs = 0.25), each = N[x])))
print(c("cor(eff, X1)" = cor(eff, X[,1]), "cor(eff, X2)" = cor(eff, X[,2])))
# We can see that eff is correlated to X2. We can confirm that the correlation is
# strongly significant.
print(c("p.value.cor(eff, X1)" = cor.test(eff, X[,1])$p.value,
        "p.value.cor(eff, X2)" = cor.test(eff, X[,2])$p.value))
# true network
G0            <- sim.network(distr)
# normalise
G0norm        <- norm.network(G0)
# Matrix GX
GX            <- peer.avg(G0norm, X)
# simulate dependent variable 
y             <- simSAR(~ -1 + eff + X + GX, Glist = G0norm, parms = c(alpha, beta), 
                        epsilon = rnorm(sum(N), sd = se))
Gy            <- y$Gy
y             <- y$y
# build dataset
dataset           <- as.data.frame(cbind(y, X, Gy, GX))
colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2")

## ----Smm7, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smmeff1 <- smmSAR(y ~ X1 + X2 || GX1 + GX2, dnetwork = distr, contextual = T,
#                       fixed.effects = T, smm.ctr  = list(R = 1, print = F),
#                       data = dataset)
# summary(out.smmeff1)

## ----Smm7a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smme1)

## ----Smm8, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smmeff2 <- smmSAR(y ~ X1 + X2 | Gy, dnetwork = distr, contextual = T,
#                       fixed.effects = T, smm.ctr  = list(R = 100, print = F),
#                       data = dataset)
# summary(out.smmeff2)

## ----Smm8a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smme2)

## ----Smm9, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smmeff3 <- smmSAR(y ~ X1 + X2, dnetwork = distr, contextual = T, fixed.effects = T,
#                       smm.ctr  = list(R = 100, print = F), data = dataset)
# summary(out.smmeff3)

## ----Smm9a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smme3)

## ----Smm10, echo = TRUE, eval = FALSE, message=FALSE--------------------------
# fMC <- function(...){
#   # Number of groups
#   M             <- 200
#   # size of each group
#   N             <- rep(30,M)
#   # individual effects
#   beta          <- c(1, 1, 1.5, 5, -3)
#   # endogenous effects
#   alpha         <- 0.4
#   # std-dev errors
#   se            <- 1
#   # network distribution
#   distr         <- runif(sum(N*(N-1)))
#   distr         <- vec.to.mat(distr, N, normalise = FALSE)
#   # covariates
#   X             <- cbind(rnorm(sum(N),0,5), rpois(sum(N),7))
#   # Groups' fixed effects
#   # We defined the groups' fixed effect as the quantile at 25% of X2 in the group
#   # This implies that the effects are correlated with X
#   eff           <- unlist(lapply(1:M, function(x)
#     rep(quantile(X[(c(0, cumsum(N))[x]+1):(cumsum(N)[x]),2], probs = 0.25), each = N[x])))
#   # true network
#   G0            <- sim.network(distr)
#   # normalise
#   G0norm        <- norm.network(G0)
#   # Matrix GX
#   GX            <- peer.avg(G0norm, X)
#   # simulate dependent variable
#   y             <- simSAR(~ -1 + eff + X + GX, Glist = G0norm,
#                           parms = c(alpha, beta), epsilon = rnorm(sum(N), sd = se))
#   Gy            <- y$Gy
#   y             <- y$y
#   # build dataset
#   dataset           <- as.data.frame(cbind(y, X, Gy, GX))
#   colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2")
#   out.smmeff1 <- smmSAR(y ~ X1 + X2 || GX1 + GX2, dnetwork = distr, contextual = T,
#                         fixed.effects = T, smm.ctr  = list(R = 1, print = F),
#                         data = dataset)
#   out.smmeff2 <- smmSAR(y ~ X1 + X2 | Gy, dnetwork = distr, contextual = T,
#                         fixed.effects = T, smm.ctr  = list(R = 100, print = F),
#                         data = dataset)
#   out.smmeff3 <- smmSAR(y ~ X1 + X2, dnetwork = distr, contextual = T, fixed.effects = T,
#                       smm.ctr  = list(R = 100, print = F), data = dataset)
#   out         <- data.frame("GX.observed"   = out.smmeff1$estimates,
#                             "Gy.observed"   = out.smmeff2$estimates,
#                             "None.observed" = out.smmeff3$estimates)
#   out
# }

## ----Smm10a, echo = TRUE, eval = FALSE, message=FALSE-------------------------
# smm.Monte.C   <- lapply(1:250, fMC)

## ----Smm10b, echo = TRUE, eval = FALSE, message=FALSE-------------------------
# Reduce('+', smm.Monte.C)/250

## ----Smm10c, echo = FALSE, eval = TRUE, message=FALSE-------------------------
print(smm$smmmc)

## ----smmsave, echo = FALSE, eval = FALSE--------------------------------------
# smm <- list(smm1 = out.smm1, smm2 = out.smm2,
#             smm3 = out.smm3, smm4 = out.smm4,
#             smme1 = out.smmeff1, smme2 = out.smmeff2,
#             smme3 = out.smmeff3, smmmc = smm$smmmc)
# save(smm, file = "~/Dropbox/Papers - In progress/Partial Network/Package/AH/PartialNetwork/datavignettes/smm.rda")

## ----Smmlogitload, echo = FALSE, eval = TRUE, message=FALSE-------------------
load(url('https://raw.githubusercontent.com/ahoundetoungan/PartialNetwork/master/datavignettes/smmlogit.rda', "rb"))

## ----smmp1, echo = TRUE, eval = FALSE-----------------------------------------
# set.seed(123)
# # Number of groups
# M        <- 100
# # size of each group
# N        <- rep(30,M)
# # covariates
# X        <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# # network formation model parameter
# rho      <- c(-0.8, 0.2, -0.1)
# # individual effects
# beta     <- c(2, 1, 1.5, 5, -3)
# # endogenous effects
# alpha    <- 0.4
# # std-dev errors
# se       <- 1
# # network
# tmp      <- c(0, cumsum(N))
# X1l      <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],1])
# X2l      <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],2])
# dist.net <- function(x, y) abs(x - y)
# X1.mat   <- lapply(1:M, function(m) {
#   matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), N[m])})
# X2.mat   <- lapply(1:M, function(m) {
#   matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), N[m])})
# Xnet     <- as.matrix(cbind("Const" = 1,
#                             "dX1"   = mat.to.vec(X1.mat),
#                             "dX2"   = mat.to.vec(X2.mat)))
# ynet     <- Xnet %*% rho
# ynet     <- c(1*((ynet + rlogis(length(ynet))) > 0))
# G0       <- vec.to.mat(ynet, N, normalise = FALSE)
# # normalise
# G0norm   <- norm.network(G0)
# # Matrix GX
# GX       <- peer.avg(G0norm, X)
# # simulate dependent variable
# y        <- simSAR(~ X +  GX, Glist = G0norm, parms = c(alpha, beta),
#                    epsilon = rnorm(sum(N), sd = se))
# Gy       <- y$Gy
# y        <- y$y
# # build dataset
# dataset           <- as.data.frame(cbind(y, X, Gy, GX))
# colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2")

## ----smmp2, echo = TRUE, eval = FALSE-----------------------------------------
# nNet      <- nrow(Xnet) # network formation model sample size
# Aobs      <- sample(1:nNet, round(0.3*nNet)) # We observed 30%
# # We can estimate rho using the gml function from the stats package
# logestim  <- glm(ynet[Aobs] ~ -1 + Xnet[Aobs,], family = binomial(link = "logit"))
# slogestim <- summary(logestim)
# rho.est   <- logestim$coefficients
# rho.var   <- slogestim$cov.unscaled # we also need the covariance of the estimator

## ----smmp3, echo = TRUE, eval = FALSE-----------------------------------------
# d.logit     <- lapply(1:M, function(x) {
#   out       <- 1/(1 + exp(-rho.est[1] - rho.est[2]*X1.mat[[x]] -
#                             rho.est[3]*X2.mat[[x]]))
#   diag(out) <- 0
#   out})

## ----Smmp4, echo = TRUE, eval = FALSE, message=FALSE--------------------------
# smm.logit   <- smmSAR(y ~ X1 + X2, dnetwork = d.logit, contextual = T,
#                       smm.ctr  = list(R = 100, print = F), data = dataset)
# summary(smm.logit)

## ----Smmp4b, echo = FALSE, eval = TRUE, message=FALSE-------------------------
print(smmlo$smm1)

## ----Smmp5, echo = TRUE, eval = FALSE, message=FALSE--------------------------
# fdist        <- function(rho.est, rho.var, M, X1.mat, X2.mat){
#   rho.est1   <- MASS::mvrnorm(mu = rho.est, Sigma = rho.var)
#   lapply(1:M, function(x) {
#   out        <- 1/(1 + exp(-rho.est1[1] - rho.est1[2]*X1.mat[[x]] -
#                              rho.est1[3]*X2.mat[[x]]))
#   diag(out)  <- 0
#   out})
# }

## ----Smmp6a, echo = TRUE, eval = FALSE, message=FALSE-------------------------
# fdist_args  <- list(rho.est = rho.est, rho.var = rho.var, M = M, X1.mat = X1.mat,
#                     X2.mat = X2.mat)
# summary(smm.logit, dnetwork = d.logit, data = dataset, .fun = fdist, .args = fdist_args,
#         sim = 500, ncores = 8) # ncores performs simulations in parallel

## ----Smmp6c, echo = FALSE, eval = TRUE, message=FALSE-------------------------
print(smmlo$smm2)

## ----smmlogitsave, echo = FALSE, eval = FALSE---------------------------------
# smm2  <- summary(smm.logit, dnetwork = d.logit, data = dataset, .fun = fdist, .args = fdist_args, sim = 500, ncores = 8)
# smmlo <- list(smm1 = smm.logit, smm2 = smm2)
# save(smmlo, file = "~/Dropbox/Papers - In progress/Partial Network/Package/AH/PartialNetwork/datavignettes/smmlogit.rda")

