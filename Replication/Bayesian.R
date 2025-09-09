library(PartialNetwork)
library(CDatanet)
library(doParallel)
library(dplyr)
library(ggplot2)
rm(list = ls())

## Work directory 
## possible root 
proot <- c("A/B/C/Working_Directory")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

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

fMC      <- function(pmis, alpha0, niter){# pmis is the proportion of non-observed pairs (i, j) in the network
  # pmis   <- 0.25
  # Matrix X
  X        <- cbind(sample(agevalue, sum(N), replace = TRUE, prob = ageprob),
                    sample(femvalue, sum(N), replace = TRUE, prob = femprob))
  
  # Matrix X centered by group
  Xc       <- do.call(rbind, lapply(1:M, function(x){
    tmp    <- colMeans(X[c(CUMN[x] + 1):CUMN[x+1],])
    t(apply(X[c(CUMN[x] + 1):CUMN[x+1],], 1, function(w) w - tmp))}))
  
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
  
  # GX
  GX       <- peer.avg(G0norm, X)
  
  # Simulate y without fixed effects
  y        <- simsar(~ X + GX, Glist = G0norm, theta = c(alpha, beta, gamma, se))
  Gy       <- y$Gy
  y        <- y$y
  
  # Observed network
  nNet      <- nrow(Xnet) # network formation model sample size
  Aobs      <- sort(sample(1:nNet, round((1 - pmis)*nNet))) #Observed entries
  AOb       <- rep(0, nNet)
  AOb[Aobs] <- 1
  AOb       <- vec.to.mat(AOb, N, normalise = FALSE)

  # Fake network
  AFa       <- sample(0:1, nNet, replace = TRUE, prob = c(0.8, 0.2))
  AFa[Aobs] <- ynet[Aobs] 
  AFa       <- vec.to.mat(AFa, N, normalise = FALSE)
  
  # logit model
  mlinks    <- list(model = "logit", 
                    mlinks.formula = ~ dX1 + dX2,
                    mlinks.data = as.data.frame(Xnet[,-1]))
  
  # Hyperparameters
  Kv         <- 2*ncol(X) + 1  # Number of exogenous explanatory variables
  hyperp     <- list("mutheta"       = rep(0,Kv),
                     "invstheta"     = diag(Kv)/100,
                     "muzeta"        = 0,
                     "invszeta"      = 1,
                     "a"             = 4.2,
                     "b"             = (4.2 - 2)*0.5)
  
  # Estimation
  dataset   <- as.data.frame(cbind(y, X)) 
  colnames(dataset) <- c("y", "X1","X2") 
  out       <- mcmcSAR(formula = y ~ X1 + X2, 
                       contextual = TRUE, 
                       G0.obs = AOb,
                       G0 = AFa, 
                       data = dataset, 
                       mlinks = mlinks, 
                       start = c(rep(0, Kv), alpha0, 1),
                       hyperparms = hyperp,
                       ctrl.mcmc = list(jumpmin = c("alpha" = 0.005, "rho" = 0.01),
                                        jumpmax = c("alpha" = 0.5, "rho" = 5),
                                        target  = c("alpha" = 0.234, "rho" = 0.234),
                                        print.level = 1),
                       iteration = niter)
}

niter  <- 1e4 # number of MCMC iterations
alph0  <- seq(-0.9, 0.9, 0.1) # Possible values for alpha0
nalph0 <- length(alph0)
fsimu  <- function(pmis, mc.cores, alph0, niter){
  cat("pmis: ", pmis, "\n")
  out  <- mclapply(alph0, function(x) fMC(pmis, x, niter), mc.cores = mc.cores)
  saveRDS(out, file = paste0("Results/Bayes.pmis=", pmis, "N=", N[1], ".RDS"))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

# 50 group of 30
M        <- 50           
N        <- rep(30,M)    
CUMN     <- c(0, cumsum(N))

fsimu(0.50, 5, alph0, niter)
fsimu(0.75, 5, alph0, niter)

# 50 group of 50
M        <- 50           
N        <- rep(50,M)    
CUMN     <- c(0, cumsum(N))

fsimu(0.50, 5, alph0, niter)
fsimu(0.75, 5, alph0, niter)


## Plot
out1 <- readRDS(file = paste0("Results/Bayes.pmis=0.5N=30.RDS"))
out2 <- readRDS(file = paste0("Results/Bayes.pmis=0.75N=30.RDS"))
out3 <- readRDS(file = paste0("Results/Bayes.pmis=0.5N=50.RDS"))
out4 <- readRDS(file = paste0("Results/Bayes.pmis=0.75N=50.RDS"))

data <- data.frame(alpha0 = factor(rep(rep(alph0, each = niter), 4)),
                   estim = unlist(c(lapply(1:nalph0, function(i) out1[[i]]$posterior$theta[,"Gy"]),
                                    lapply(1:nalph0, function(i) out2[[i]]$posterior$theta[,"Gy"]),
                                    lapply(1:nalph0, function(i) out3[[i]]$posterior$theta[,"Gy"]),
                                    lapply(1:nalph0, function(i) out4[[i]]$posterior$theta[,"Gy"]))),
                   simu  = rep(1:niter, 4*nalph0),
                   model = rep(1:4, each = niter*nalph0)) %>%
  mutate(model = factor(model, labels = c(
    "paste('Share of missing pairs = 50%, ', M, '= 50, ', N[m], '= 30')",
    "paste('Share of missing pairs = 75%, ', M, '= 50, ', N[m], '= 30')",
    "paste('Share of missing pairs = 50%, ', M, '= 50, ', N[m], '= 50')",
    "paste('Share of missing pairs = 75%, ', M, '= 50, ', N[m], '= 50')"
  )))

(graph <- ggplot(data = data %>% filter(simu <= 5e3), aes(x = simu, y = estim, colour = alpha0)) + 
    geom_line() +
    xlab("MCMC Iteration") +
    ylab("Peer Effect Estimate") +
    facet_wrap(~ model, labeller = label_parsed) + 
    theme_bw(base_size = 12, base_family = "Palatino")+
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(-0.9, 0.9), breaks = round(seq(-0.9, 0.9, 0.3), 1))) 

ggsave("Bayesian.pdf", plot = graph, device = "pdf", width = 8, height = 5)
