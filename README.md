# An R package for estimating Peer Effects Using Partial Network Data

The **PartialNetwork** package includes all functions for the replication of the results in Boucher and Houndetoungan (2023). The exact replication codes are located in the folder [**test**](https://github.com/ahoundetoungan/PartialNetwork/tree/master/test). Below, we also provide detailed examples on how to use the estimators described in the paper.

## Installation
### CRAN version
**PartialNetwork** can be directly installed from CRAN.
```R
install.packages("PartialNetwork")
```

### GitHub version
It may be possible that we updated the package without submitting the new version to CRAN. The latest version (*but not necessary stable*) of **PartialNetwork** can be installed from this GitHub repos. 
```R
remotes::install_github("ahoundetoungan/PartialNetwork", build_vignettes = TRUE)
```
## How to use the `PartialNetwork` package?
See our [vignette in pdf](https://nbviewer.jupyter.org/github/ahoundetoungan/PartialNetwork/blob/master/doc/PartialNetwork_vignette.pdf).
<!--## Instrumental Variable procedure

We provide the function `sim.IV(dnetwork, X, y, replication, power)` where `dnetwork` is the network linking probabilities, `X` is a matrix of covariates, `y` (optional) is the vector of outcome, `replication` (optional, default = 1) is the number of replication, and `power` (optional, default = 1) is the number of powers of the interaction matrix used to generate the instruments. The function outputs a proxy for Gy and simulated instruments. See the help file (`? sim.IV`) of the function for details. The following code provides an example using a single group of network.

```R
# initialize parameters
# number of groups
M             <- 5
# size of the group
N             <- 100      
# value of lambda (precision parameter for the network formation model)
lambda        <- 1   
# individual effects
beta          <- c(2, 1, 1.5) 
# contextual effects
gamma         <- c(5, -3) 
# endogenous effect
alpha         <- 0.4
# std-dev errors
se            <- 1

# network probabilities
Prob          <- list()

for (m in 1:M) {
  # heterogeneity of the linking probabilities
  c           <- rnorm(N*N, 0, 1) 
  # generate linking probabilities
  Pm          <- matrix(exp(c / lambda) / (1 + exp(c / lambda)), N)
  # no self-link
  diag(Pm)    <- 0 
  Prob[[m]]   <- Pm
}

# generate data
X             <- matrix(data =  NA, nrow = 0, ncol = 2)
y             <- c()

for (m in 1:M) {
  # generate the 'observed network'
  Gm          <- sim.network(Prob[[m]]) 
  rs          <- rowSums(Gm)
  rs[rs == 0] <- 1
  # row-normalize
  Wm          <- Gm/rs
  # covariates
  Xm          <- cbind(rnorm(N,0,5),rpois(N,6)) 
  # endogenous variable, no contextual effect
  ym          <- solve(diag(N) - alpha * Wm) %*% (cbind(rep(1, N), Xm) %*% beta + rnorm(N,0,se)) 
  y           <- c(y, ym)
  X           <- rbind(X, Xm)
}

# generate instruments 
instr         <- sim.IV(Prob, X, y, replication = 1, power = 2)

GY1c1         <- instr[[1]]$G1y       # proxy for Gy (draw 1)
GXc1          <- instr[[1]]$G1X[,,1]  # proxy for GX (draw 1)
G2Xc1         <- instr[[1]]$G1X[,,2]  # proxy for GGX (draw 1)
GXc2          <- instr[[1]]$G2X[,,1]  # proxy for GX (draw 2)
G2Xc2         <- instr[[1]]$G2X[,,2]  # proxy for GGX (draw 2)
```
Once the instruments are generated, the estimation can be performed using standard tools, e.g. the function `ivreg` from the [**AER**](https://cran.r-project.org/package=AER) package. For example:
  ```R
# build dataset
# keep only instrument constructed using a different draw than the one used to proxy Gy
dataset           <- as.data.frame(cbind(y,X,GY1c1,GXc2,G2Xc2)) 
# rename variables
colnames(dataset) <- c("y","X1","X2","Gy1","Z1","Z2","ZZ1","ZZ2") 
library(AER)
results           <- ivreg(y ~ X1 + X2 + Gy1 | X1 + X2 + Z1 + Z2 + ZZ1 + ZZ2, data = dataset)
summary(results)
```

## Bayesian estimator

The Bayesian estimator is neatly packed in the function `mcmcSAR(formula, contextual, start, G0, hyperparms, iteration, ctrl.mcmc, data)`, where `formula` is the model equation, `contextual` indicates if the model has contextual effects, `start` (optional) is the parameter initialization, `G0` (optional) is the starting of the network,  `hyperparam` specify the prior distributions, including the network linking probabilities, `iterations` (optional, default = 2000) is the number of MCMC steps to be performed, `ctrl.mcmc` set some controls for the MCMC and `data` contains the data (if not specified R will search the variables in the global environment). See the help (`? mcmcSAR`) file for a complete description. Below, we provide a simple example using simulated data.

### Simulate data
```R
# Number of groups
M             <- 100
# size of each group
N             <- rep(50,M)
# precision parameter for the network formation process
lambda        <- 1 

G             <- list()

# individual effects
beta          <- c(2,1,1.5) 
# contextual effects
gamma         <- c(5,-3) 
# endogenous effect
alpha         <- 0.4
# std-dev errors
se            <- 1 

prior         <-list()

## generate network probabilities
for (m in 1:M) {
  Nm          <- N[m]
  c           <- rnorm(Nm*Nm,0,1)
  # linking probabilities
  Prob        <- matrix(exp(c/lambda)/(1+exp(c/lambda)),Nm) 
  # no self-link
  diag(Prob)  <- 0 
  prior[[m]]  <-Prob
}

## generate data
# covariates
X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# dependent variable
y             <- c()

for (m in 1:M) {
  Nm          <- N[m]
  # true network
  Gm          <- matrix(runif(Nm^2),Nm,Nm) < prior[[m]] 
  # no self-link
  diag(Gm)    <- rep(0,Nm) 
  G[[m]]      <- Gm
  rsm         <- rowSums(Gm)
  rsm[rsm==0] <- 1
  # normalize
  Gm          <- Gm/rsm 
  # rows index of group m
  r2          <- sum(N[1:m])
  r1          <- r2 - Nm + 1
  # contextual effect
  Xm          <- X[r1:r2,]
  GXm         <- Gm %*% Xm
  y[r1:r2]    <- solve(diag(Nm)-alpha*Gm) %*% (cbind(rep(1,Nm),Xm) %*% beta + GXm %*% gamma + rnorm(Nm,0,se)) 
}
```

### Estimate the model on simulated data
```R
# set the hyperparameter
# the hyperparameter is a list
hyperparms    <- list("dnetwork" = prior) 

# launch the MCMC
ctrl          <- list(print.level = 2)
out           <- mcmcSAR(y ~ X | X, hyperparms = hyperparms, ctrl.mcmc = ctrl)

# the packahe has summary and plot methods of class S3
# summary
summary(out)
# plot simulations
plot(out)
# plot density functions
plot(out, plot.type = "dens")
```
## ARD, Breza et al. (2020)

### Simulation procedure
The data is simulated following a procedure similar to the one in Breza et al. (2020). One notable exception is how we attribute traits to individuals. The code uses the functions `rvMF` (random variable generator von Mises-Fisher (vMF)), `dvMF` (density vMF), `CpvMF` (Normalizing constant, vMF), `sim.dnetwork` (linking probabilities) and `sim.network` (drawn network from probabilities). See the associated help files for details. Below, we provide an example using only one group.

```R
# Sample size
N       <- 500 

# ARD parameters
genzeta <- 1
mu      <- -1.35
sigma   <- 0.37
K       <- 12    # number of traits
P       <- 3     # Sphere dimension 


# Generate z (spherical coordinates)
genz    <- rvMF(N,rep(0,P))

# Generate nu  from a Normal distribution with parameters mu and sigma (The gregariousness)
gennu   <- rnorm(N,mu,sigma)

# compute degrees
gend    <- N*exp(gennu)*exp(mu+0.5*sigma^2)*exp(logCpvMF(P,0) - logCpvMF(P,genzeta))

# Link probabilities
Probabilities <- sim.dnetwork(gennu,gend,genzeta,genz) 

# Adjacency matrix
G        <- sim.network(Probabilities)

# Generate vk, the trait location
genv     <- rvMF(K,rep(0,P))

# set fixed some vk  distant
genv[1,] <- c(1,0,0)
genv[2,] <- c(0,1,0)
genv[3,] <- c(0,0,1)

# eta, the intensity parameter
geneta   <-abs(rnorm(K,2,1))

# Build traits matrix
densityatz       <- matrix(0,N,K)
for(k in 1:K){
  densityatz[,k] <- dvMF(genz,genv[k,]*geneta[k])
}
trait       <- matrix(0,N,K)

for(k in 1:K){
  trait[,k] <- densityatz[,k]>sort(densityatz[,k],decreasing = T)[runif(1,0.05*N,0.25*N)]
}
# print a percentage of people having a trait
colSums(trait)*100/N
  
# Build ADR
ARD         <- G %*% trait
  
# generate b
genb        <- numeric(K)
for(k in 1:K){
   genb[k]  <- sum(G[,trait[,k]==1])/sum(G)
}

```
### Estimate the model on simulated data

We present a simple function wrapping, `mcmcARD`, for the estimation procedure proposed by Breza et al. (2020). For specific information on the function, see the help file.

```R
# initialization 
d0     <- exp(rnorm(N)); b0 <- exp(rnorm(K)); eta0 <- rep(1,K);
zeta0  <- 05; z0 <- matrix(rvMF(N,rep(0,P)),N); v0 <- matrix(rvMF(K,rep(0,P)),K)
  
# We should fix some vk and bk
vfixcolumn      <- 1:5
bfixcolumn      <- c(3, 5)
b0[bfixcolumn]  <- genb[bfixcolumn]
v0[vfixcolumn,] <- genv[vfixcolumn,]
start           <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0, "zeta" = zeta0)
  
# MCMC
out   <- mcmcARD(Y = ARD, traitARD = trait, start = start, fixv = vfixcolumn, consb = bfixcolumn, iteration = 5000)

# fit network distribution using the estimates from mcmcARD
dist  <- fit.dnetwork(out, traitARD = trait)
```-->
