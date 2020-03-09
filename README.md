# PartialNetwork
Estimating Peer Effects Using Partial Network Data

## How to install
```R
ddevtools::install_github("PartialNetwork")
```

## Instrumental Variable procedure

We provide the function *instrument(distr,X,Y,S,pow)* where *distr* is the network linking probabilities, *X* is a matrix of covariates, *Y* (optional) is the vector of outcome, *S* (optional, default=2) is the number of network draws used, and *pow* is the number of powers of the interaction matrix used to generate the instruments. The function outputs a proxy for Gy and simulated instruments. See the help file of the function for details. The following code provides an example using a single group.

```{c}
# initialize parameters
N <- 50
lambda <- 1
beta <- c(2, 1, 1.5)
gamma <- c(5,-3)
alpha <- 0.4
se <- 1
c <- rnorm(N*N,0,1)

# network probabilities
Probabilities <- matrix(exp(c / lambda) / (1 + exp(c / lambda)), N)
diag(Probabilities) <- 0

# generate data
G<-Graph(Probabilities)
rs<-rowSums(G)
rs[rs==0]<-1
W<-G/rs
X <- cbind(rnorm(N,0,5),rpois(N,6))
Y <- solve(diag(rep(1,N))-alpha*W)%*%(cbind(rep(1,N), X)%*%beta + rnorm(N,0,se))

# generate instruments
instr1 <- instruments(Probabilities, X, Y, S=2, pow=2)
GY1c <- instr1$GY  # proxy for Gy (draw 1)
GXc0 <-instr1$GX[[1]][,,1] # proxy for GX (draw 1)
G2Xc0 <-instr1$GX[[1]][,,2]  # proxy for GGX (draw 1)
GXc <- instr1$GX[[2]][,,1]  # proxy for GX (draw 2)
G2Xc <- instr1$GX[[2]][,,2]  # proxy for GGX (draw 2)
```

## Bayesian estimator

The estimator is neatly packed in the function *peerMCMC(y,X,para,hyperparams,iteractions)*. See the help file for a complete description. Below, we provide a simple example using simulated data.

### Simulate data
```{r}
M <- 100 #Number of groups
N <- rep(50,M)
lambda <- 1

bayesian.estimate <- matrix(0,0,7)
colnames(bayesian.estimate) <- c("Intercept", paste0("X",1:2), paste0("GX",1:2), "alpha", "sigma^2")
MCMC.iteration <- 2000

G<-list()
X<-list()
y<-list()

beta <- c(2,1,1.5)
gamma <- c(5,-3)
alpha <- 0.4
se <- 1

out   <- list()
prior <-list()
  
## generate network probabilities
for (m in 1:M) {
  xx <- sample(1:N[m],3*N[m],TRUE)
  yy <- sample(1:N[m],3*N[m],TRUE)
  c <- rnorm(N[m]*N[m],0,1)
  Probabilities <- matrix(exp(c/lambda)/(1+exp(c/lambda)),N[m])
  diag(Probabilities) <- 0
  prior[[m]]<-Probabilities
}
  
## generate data
for (m in 1:M) {
  Gm <- matrix(runif(N[m]^2),N[m],N[m]) < prior[[m]]
  diag(Gm) <- rep(0,N[m])
  G[[m]] <- Gm
  rsm <- rowSums(Gm)
  rsm[rsm==0] = 1
  Gm = Gm/rsm
  Xm <- cbind(rnorm(N[m],0,5),rpois(N[m],7))
  if (m %% 2 == 0) Xm[,1] <- 0
  if (m %% 2 == 1) Xm[,2] <- 7
  X[[m]] <- Xm
  GXm <- Gm%*%Xm
  y[[m]] <- solve(diag(rep(1,N[m]))-alpha*Gm)%*%(cbind(rep(1,N[m]),Xm)%*%beta+GXm%*%gamma+rnorm(N[m],0,se))
}
```
### Estimate the model on simulated data
```{r}
Kv <- 2*ncol(X[[1]]) + 1 # number of parameters

## set the model's options
hyperparms <- list("network" = prior, # network prior distribution
                     "theta0" = rep(0,Kv), # initial parameter value
                     "invsigmatheta" = diag(Kv)/100, # variance-covariance matrix
                     "zeta0" = 0, # ???
                     "insigma2zeta" = 2, # ???
                     "a" = 4.2, #hyperparameter for sigma2
                     "b" = (4.2 - 2))  #hyperparameter for sigma2
  
### launch the MCMC
out <- peerMCMC(y, X, c(beta,gamma,0.4,se), hyperparms, iteration = MCMC.iteration)
```

## ARD, Breza et al. (2020)

### Simulation procedure
The data is simulated following a procedure similar to the one in Breza et al. (2020). One notable exception is how we attribute traits to individuals. The code uses the functions *rvMF*, *dvMF*, *Prob* and *Graph*. See the associated help files for details.

```{r}

N <- 250
genzeta <- 1.5
mu <- -1.25
sigma <- 0.37
K <- 12
P <- 3

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

We present a simple function wrapping, *updateGP*, for the estimation procedure proposed by Breza et al. (2020). For specific information on the function, see the help file.

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
