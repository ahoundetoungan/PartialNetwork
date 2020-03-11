# PartialNetwork
Estimating Peer Effects Using Partial Network Data

This package includes all functions for the replication of the results in Boucher and Houndetoungan (2020). The exact replication codes are located in (DIRECTORY). Below, we also provide detailed examples on how to use the estimators described in the paper.

## How to install
```R
ddevtools::install_github("PartialNetwork")
```

## Instrumental Variable procedure

We provide the function `instrument(distr,X,Y,S,pow)` where `distr` is the network linking probabilities, `X` is a matrix of covariates, `Y` (optional) is the vector of outcome, `S` (optional, default=2) is the number of network draws used, and `pow` (optional) is the number of powers of the interaction matrix used to generate the instruments. The function outputs a proxy for Gy and simulated instruments. See the help file of the function for details. The following code provides an example using a single group.

```{c}
# initialize parameters
N <- 50 # size of the group
lambda <- 1 # value of lambda (precision parameter for the network formation model)
beta <- c(2, 1, 1.5) # individual effects
gamma <- c(5,-3) # contextual effects
alpha <- 0.4 # endogenous effect
se <- 1 # std-dev errors
c <- rnorm(N*N,0,1) # heterogeneity of the linking probabilities

# network probabilities
Probabilities <- matrix(exp(c / lambda) / (1 + exp(c / lambda)), N) # generate linking probabilities
diag(Probabilities) <- 0 # no self-link

# generate data
G<-Graph(Probabilities) # generate the 'observed network'
rs<-rowSums(G)
rs[rs==0]<-1
W<-G/rs # row-normalize
X <- cbind(rnorm(N,0,5),rpois(N,6)) # covariates
Y <- solve(diag(rep(1,N))-alpha*W)%*%(cbind(rep(1,N), X)%*%beta + rnorm(N,0,se)) # endogenous variable, no contextual effect

# generate instruments ### ARISTIDE ???? Could you check this is correct?
instr1 <- instruments(Probabilities, X, Y, S=2, pow=2)
GY1c1 <- instr1$GY  # proxy for Gy (draw 1)
GXc1 <-instr1$GX[[1]][,,1] # proxy for GX (draw 1)
G2Xc1 <-instr1$GX[[1]][,,2]  # proxy for GGX (draw 1)
GXc2 <- instr1$GX[[2]][,,1]  # proxy for GX (draw 2)
G2Xc2 <- instr1$GX[[2]][,,2]  # proxy for GGX (draw 2)
```
Once the instruments are generated, the estimation can be performed using standard tools, e.g. the function `ivreg` from the AER package (required by PartialNetwork). For example:
```{r}
dataset <- as.data.frame(cbind(Y,X,GY1c1,GX2,G2X2)) # build dataset, keep only instrument constructed using a different draw than than the one used to proxy Gy
colnames(dataset) <- c("Y","X1","X2","Gy1","Z1","Z2","ZZ1","ZZ2") # rename variables
results <- ivreg(Y~ X1 + X2 + Gy1 | X1 + X2 + Z1 + Z2 + ZZ1 + ZZ2, data=dataset)
```

## Bayesian estimator

Our Bayesian estimator is neatly packed in the function `peerMCMC(y,X,para,hyperparams,iteractions)`, where `y` is the list (lenght=M) of endogenous variables, `X` is the list of covariates, `hyperparam` specify the prior distributions, including the network linking probabilities, and `iterations` is the number of MCMC steps to be performed. See the help file for a complete description. Below, we provide a simple example using simulated data.

### Simulate data
```{r}

## initialize
M <- 100 #Number of groups
N <- rep(50,M) # size of each group
lambda <- 1 # precision parameter for the network formation process
bayesian.estimate <- matrix(0,0,7)
colnames(bayesian.estimate) <- c("Intercept", paste0("X",1:2), paste0("GX",1:2), "alpha", "sigma^2")
MCMC.iteration <- 2000

G<-list()
X<-list()
y<-list()

beta <- c(2,1,1.5) # individual effects
gamma <- c(5,-3) # contextual effects
alpha <- 0.4 # endogenous effect
se <- 1 # std-dev errors

out   <- list()
prior <-list()
  
## generate network probabilities
for (m in 1:M) {
  c <- rnorm(N[m]*N[m],0,1)
  Probabilities <- matrix(exp(c/lambda)/(1+exp(c/lambda)),N[m]) # linking probabilities
  diag(Probabilities) <- 0 # no self-link
  prior[[m]]<-Probabilities
}
  
## generate data
for (m in 1:M) {
  Gm <- matrix(runif(N[m]^2),N[m],N[m]) < prior[[m]] # observed network
  diag(Gm) <- rep(0,N[m]) # no self-link
  G[[m]] <- Gm
  rsm <- rowSums(Gm)
  rsm[rsm==0] = 1
  Gm = Gm/rsm # normalize
  Xm <- cbind(rnorm(N[m],0,5),rpois(N[m],7)) # covariates
  if (m %% 2 == 0) Xm[,1] <- 0  ## ARISTIDE ???
  if (m %% 2 == 1) Xm[,2] <- 7  ## ARISTIDE ???
  X[[m]] <- Xm
  GXm <- Gm%*%Xm # contextual effect
  y[[m]] <- solve(diag(rep(1,N[m]))-alpha*Gm)%*%(cbind(rep(1,N[m]),Xm)%*%beta+GXm%*%gamma+rnorm(N[m],0,se)) # endogenous variable
}
```

### Estimate the model on simulated data
```{r}
Kv <- 2*ncol(X[[1]]) + 1 # number of parameters

## set the model's options
hyperparms <- list("network" = prior, # network prior distribution
                     "theta0" = rep(0,Kv), # ARISTIDE ???
                     "invsigmatheta" = diag(Kv)/100, # ARISTIDE ???
                     "zeta0" = 0, # ARISTIDE ???
                     "insigma2zeta" = 2, # ARISTIDE ???
                     "a" = 4.2, #ARISTIDE ???
                     "b" = (4.2 - 2))  #ARISTIDE ???
  
### launch the MCMC
out <- peerMCMC(y, X, c(beta,gamma,0.4,se), hyperparms, iteration = MCMC.iteration)
```

## ARD, Breza et al. (2020)

### Simulation procedure
The data is simulated following a procedure similar to the one in Breza et al. (2020). One notable exception is how we attribute traits to individuals. The code uses the functions `rvMF` (random variable generator vonMises-Fisher (vMF)), `dvMF` (density vMF), `CpvMF` (Normalizing constant, vMF), `Prob` (linking probabilities) and `Graph` (drawn network from probabilities). See the associated help files for details. Below, we provide an example using only one group.

```{r}

N <- 250 # size of the group
genzeta <- 1.5
mu <- -1.25
sigma <- 0.37
K <- 12 # number of traits
P <- 3 # number of dimensions of the hypersphere

#1- Generate z
genz <- rvMF(N, rep(0, P))

#2- Genetate nu  from a Normal distribution with parameters mu and sigma
gennu <- rnorm(N, mu, sigma)

#3- Generate a graph G

#compute d
gend <- N * exp(gennu) * exp(mu + 0.5 * sigma ^ 2) * (CpvMF(P, 0) / CpvMF(P, genzeta))

#Network
Probabilities <- Prob(gennu, gend, genzeta, genz) 
G <- Graph(Probabilities)

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
densityatz <- matrix(0, N, K)
for (k in 1:K) {
  densityatz[, k] <- dvMF(genz, genv[k, ] * geneta[k])
}
    
## generate traits
trait <- matrix(0, N, K)
NK <- floor(runif(K,0.8,0.95)*colSums(densityatz)/unlist(lapply(1:K, function(w){max(densityatz[,w])})))
    
for (k in 1:K) {
  trait[,k]<-rbinom(N,1,NK[k]*densityatz[,k]/sum(densityatz[,k]))
}

#5 contruct ADR
ARD <- G %*% trait

#generate b
genb <- numeric(K)
for (k in 1:K) {
  genb[k] <- sum(G[, trait[, k] == 1]) / sum(G)
}

```
### Estimate the model on simulated data

We present a simple function wrapping, `updateGP`, for the estimation procedure proposed by Breza et al. (2020). For specific information on the function, see the help file.

```{r}
# initionalization
d0 <- exp(rnorm(N))
b0 <- exp(rnorm(K))
eta0 <- rep(1, K)
    
zeta0 <- 05
z0 <- matrix(rvMF(N, rep(0, P)), N)
v0 <- matrix(rvMF(K, rep(0, P)), K)

#We should fix one bk for identification
vfixcolumn <- 1:5
bfixcolumn <- 3
b0[bfixcolumn] <- genb[bfixcolumn]
v0[vfixcolumn, ] <- genv[vfixcolumn, ]
    
#Initialization
Begin<-list(z0,v0,d0,b0,eta0,zeta0)
hparms = c(0, 1, 0, 1, 5, 0.5, 1, 1)
Iter = 5000

#Update the parameters
Gparmsupdate <- updateGP(ARD, trait, Begin, vfixcolumn, bfixcolumn, Iter)
```
