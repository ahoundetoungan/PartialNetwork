library(PartialNetwork)
library(CDatanet)
library(doParallel)
rm(list = ls())
## Work directory 
## possible root 
proot <- c("~/Dropbox/Academy/1.Papers/Partial Network/Simulations/Monte Carlo",
           "~/PartialNetwork/Monte Carlo")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)
Rcpp::sourceCpp("CppFunctions.cpp")

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
rho      <- c(-2.346986, 0.404115, -0.699205)
se       <- 0.707

fMC      <- function(pmisclass){# pmisclass includes the proportion of false positive and negative
  type   <- ifelse(all(c("p0", "p1") %in% names(pmisclass)), 3, ifelse(names(pmisclass) == "p1", 1, 2))
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
  
  # Observed network
  nNet      <- nrow(Xnet) # network formation model sample size
  ynF       <- ifelse(runif(nNet) < ifelse(is.na(pmisclass["p1"]), 0, pmisclass["p1"]) & ynet == 1, 0, 
                      ifelse(runif(nNet) <  ifelse(is.na(pmisclass["p0"]), 0, pmisclass["p0"]) & ynet == 0, 1, ynet))# The fake ynet
  AOb       <- vec.to.mat(ynF, N, normalise = FALSE)
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
  
  # Estimation the peer effect model using the observed network 
  # Without fixed effects
  W           <- solve(crossprod(cbind(1, X, GObX, GGObX))/sum(N))
  gmm         <- smmSAR(y ~ X1 + X2|GOby|GObX1 + GObX2, contextual = T, dnetwork = AOb, cond.var = FALSE,
                        W = W, smm.ctr  = list(R = 1, print = F), data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmm["Gy"]   <- ifelse(gmm["Gy"] < -1, -1, ifelse(gmm["Gy"] > 1, 1, gmm["Gy"]))
  names(gmm)  <- paste0("gmm", names(gmm))
  
  # With fixed effects
  Wf          <- solve(crossprod(cbind(Xc, GObXc, GGObXc))/sum(N))
  gmmf        <- smmSAR(yf ~ X1 + X2|GObyf|GObX1 + GObX2, contextual = T, dnetwork = AOb, cond.var = FALSE,
                        fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 1, print = F), 
                        data = dataset)$estimates
  # This estimator has a closed-form expression and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  gmmf["Gy"]  <- ifelse(gmmf["Gy"] < -1, -1, ifelse(gmmf["Gy"] > 1, 1, gmmf["Gy"]))
  names(gmmf) <- paste0("gmmf", names(gmmf))
  
  # Estimation of the network distribution
  estLog    <- flogistmisclassify(c(0.999 - pmisclass, rho), y = ynF, X = Xnet, type = type)
  p00.est   <- ifelse(type == 1, 1, estLog$estimate[1]) # Specificity estimate
  p11.est   <- ifelse(type == 1, estLog$estimate[1], ifelse(type == 2, 1, estLog$estimate[2])) # Sensitivity estimate
  rho.est   <- estLog$estimate[-c(1, ifelse(type == 3, 2, 1))]    
  hpl       <- c(1/(1 + exp(-as.matrix(Xnet)%*%rho.est)))
  d.logit   <- ifelse(ynF == 0, (1 - p11.est)*hpl/((1 - p11.est)*hpl + p00.est*(1 - hpl)),
                      p11.est*hpl/(p11.est*hpl + (1 - p00.est)*(1 - hpl)))
  d.logit   <- vec.to.mat(d.logit, N, normalise = FALSE)

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
  # These estimators have closed-form expressions and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  smm1["Gy"]   <- ifelse(smm1["Gy"] < -1, -1, ifelse(smm1["Gy"] > 1, 1, smm1["Gy"]))
  smm1f["Gy"]  <- ifelse(smm1f["Gy"] < -1, -1, ifelse(smm1f["Gy"] > 1, 1, smm1f["Gy"]))
  names(smm1)  <- paste0("smm1", names(smm1))
  names(smm1f) <- paste0("smm1f", names(smm1f))
  
  # SMM with Gy observed and GX unobserved
  smm2       <- smmSAR(y ~ X1 + X2|Gy, contextual = T, dnetwork = d.logit, 
                       W = W, smm.ctr  = list(R = 100, print = F), data = dataset)$estimates
  smm2f      <- smmSAR(yf ~ X1 + X2|Gyf, contextual = T, dnetwork = d.logit, 
                       fixed.effects = TRUE, W = Wf, smm.ctr  = list(R = 100, print = F), 
                       data = dataset)$estimates
  # These estimators have closed-form expressions and may be lower than -1 or higher than 1
  # We constrain it to be coherent with the estimators that do not have a closed-form expression.
  smm2["Gy"]   <- ifelse(smm2["Gy"] < -1, -1, ifelse(smm2["Gy"] > 1, 1, smm2["Gy"]))
  smm2f["Gy"]  <- ifelse(smm2f["Gy"] < -1, -1, ifelse(smm2f["Gy"] > 1, 1, smm2f["Gy"]))
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
  
  c(gmm, gmmf, smm1, smm1f, smm2, smm2f, smm3, smm3f, smm4, smm4f)
}

fsimu  <- function(pmisclass, mc.cores){
  out  <- do.call(rbind, mclapply(1:1e3, function(i) fMC(pmisclass), mc.cores = mc.cores))
  saveRDS(out, file = paste0("Results/logit.misclass.p0=", ifelse(is.na(pmisclass["p0"]), 0, pmisclass["p0"]), 
                             ".p1=", ifelse(is.na(pmisclass["p1"]), 0, pmisclass["p1"]), ".RDS"))
}
set.seed(123)
fsimu(c("p0" = 0.05), 15)
fsimu(c("p0" = 0.10), 15)
fsimu(c("p0" = 0.25), 15)
fsimu(c("p0" = 0.50), 15)
fsimu(c("p1" = 0.05), 15)
fsimu(c("p1" = 0.10), 15)
fsimu(c("p1" = 0.25), 15)
fsimu(c("p1" = 0.50), 15)
fsimu(c("p0" = 0.05, "p1" = 0.05), 15)
fsimu(c("p0" = 0.10, "p1" = 0.10), 15)
fsimu(c("p0" = 0.25, "p1" = 0.25), 15)
fsimu(c("p0" = 0.5, "p1" = 0.5), 15)

out <- cbind(t(apply(readRDS("Results/logit.misclass.p0=0.05.p1=0.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.1.p1=0.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.25.p1=0.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.5.p1=0.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.p1=0.05.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.p1=0.1.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.p1=0.25.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.p1=0.5.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.05.p1=0.05.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.1.p1=0.1.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.25.p1=0.25.RDS"), 2, function(x) c(mean(x), sd(x)))),
             t(apply(readRDS("Results/logit.misclass.p0=0.5.p1=0.5.RDS"), 2, function(x) c(mean(x), sd(x)))))

colnames(out) <- paste0(rep(c("mean ", "sd "), ncol(out)/2),
                        rep(c("p0=5% p1=0%", "p0=10% p1=0%", "p0=25% p1=0%", "p0=50% p1=0%", 
                              "p0=0% p1=5%", "p0=0% p1=10%", "p0=0% p1=25%", "p0=0% p1=50%",
                              "p0=5% p1=5%", "p0=10% p1=10%", "p0=25% p1=25%", "p0=50% p1=50%"), each = 2))

write.csv(out, file = "Results/logit.misclassified.csv")

# Graph
library(ggplot2)
library(dplyr)
FPR  <- c(0.05, 0.1, 0.25, 0.5, rep(0, 4), 0.05, 0.1, 0.25, 0.5)
FNR  <- c(rep(0, 4), 0.05, 0.1, 0.25, 0.5, 0.05, 0.1, 0.25, 0.5)
FPRL <- c("5%", "10%", "25%", "50%", rep("0%", 4), "5%", "10%", "25%", "50%")
FNRL <- c(rep("0%", 4), "5%", "10%", "25%", "50%", "5%", "10%", "25%", "50%")
data <- do.call(rbind, lapply(1:length(FPR), function(k) {
  est  <- apply(readRDS(paste0("Results/logit.misclass.p0=", FPR[k], ".p1=", FNR[k], ".RDS")),
                2, function(x) c(mean(x), quantile(x, prob = c(0.025, 0.975))))
  data.frame(FPR   = FPRL[k],
             FNR   = FNRL[k],
             MClas = k,
             spec  = rep(0:4, each = 2),
             model = c(rep("Classical IV", 2),
                       rep("SGMM: Gy, GX observed", 2),
                       rep("SGMM: Gy observed, GX unobserved", 2),
                       rep("SGMM: Gy unobserved, GX observed", 2),
                       rep("SGMM: Gy, GX unobserved", 2)),
             FE    = rep(c(FALSE, TRUE),  5),
             coef  = est[1, c("gmmGy", "gmmfGy", "smm1Gy", "smm1fGy", 
                              "smm2Gy", "smm2fGy", "smm3Gy", "smm3fGy", 
                              "smm4Gy", "smm4fGy")],
             IC1   = est[2, c("gmmGy", "gmmfGy", "smm1Gy", "smm1fGy", 
                              "smm2Gy", "smm2fGy", "smm3Gy", "smm3fGy", 
                              "smm4Gy", "smm4fGy")],
             IC2   = est[3, c("gmmGy", "gmmfGy", "smm1Gy", "smm1fGy", 
                              "smm2Gy", "smm2fGy", "smm3Gy", "smm3fGy", 
                              "smm4Gy", "smm4fGy")])
}))

data   <- data %>% mutate(Model = factor(spec, labels = unique(model)),
                          MC    = paste0(FPR, "\n", FNR),
                          MC    = factor(MC, levels = unique(MC[order(MClas)])))

(NOFE  <- ggplot(data %>% filter(MClas %in% c(3, 4, 7, 8, 11, 12),
                                 FE == FALSE, spec %in% c(0:4)), aes(x = MC, colour = Model)) +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "gray") +
    geom_errorbar(width=.2, aes(ymin = IC1, ymax = IC2),
                  position = position_dodge(width = 0.3)) +
    geom_point(aes(y = coef, shape = Model), position = position_dodge(width = 0.3)) +
    theme_bw() +
    annotate("text", x = 0.2, y = alpha, label = bquote(alpha[0] == .(round(alpha, 3))), 
             hjust = -0.1, vjust = -0.2, size = 3, color = "#333") +
    xlab("False positive rate (first row) and false negative rate (second row)") + ylab("Peer effect estimate") +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)))
ggsave("mc_mclasNOFE.pdf", plot = NOFE, device = "pdf", width = 9, height = 4)

(WIFE  <- ggplot(data %>% filter(MClas %in% c(3, 4, 7, 8, 11, 12),
                                 FE == TRUE, spec %in% c(0:4)), aes(x = MC, colour = Model)) +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "gray") +
    geom_errorbar(width=.2, aes(ymin = IC1, ymax = IC2),
                  position = position_dodge(width = 0.3)) +
    geom_point(aes(y = coef, shape = Model), position = position_dodge(width = 0.3)) +
    theme_bw() +
    annotate("text", x = 0.2, y = alpha, label = bquote(alpha[0] == .(round(alpha, 3))), 
             hjust = -0.1, vjust = -0.2, size = 3, color = "#333") +
    xlab("False positive rate (first row) and false negative rate (second row)") + ylab("Peer effect estimate") +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)))
ggsave("mc_mclasWIFE.pdf", plot = WIFE, device = "pdf", width = 9, height = 4)

