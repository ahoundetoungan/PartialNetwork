#####################################Headline#####################################
rm(list = ls())
setwd("~/Dropbox/ARD")
setwd("~/ARD")
library(Rcpp)                         #Install Rcpp if not already done. (also the libraries RcppArmadillo, RcppEigen and BH)
sourceCpp("Simulations/CppFiles/GX.not.observed2.cpp")
library(ggplot2)
##################################################################################

####################################Simulations###################################
set.seed(123)

M <- 100 #Number of groups
N <- rep(50,M)
lambda <- 1

bayesian.estimate <- matrix(0,0,7)
colnames(bayesian.estimate) <- c("Intercept", paste0("X",1:2), paste0("GX",1:2), "alpha", "sigma^2")

MC.iteration <- 2
MCMC.iteration <- 2000

G<-list()
X<-list()
y<-list()

beta <- c(2,1,1.5)
gamma <- c(5,-3)
alpha <- 0.4
se <- 1

#Iteration <- 1000

out   <- list()
prior <-list()
t1 <- Sys.time()
for (l in 1:MC.iteration) {
  for (m in 1:M) {
    xx <- sample(1:N[m],3*N[m],TRUE)
    yy <- sample(1:N[m],3*N[m],TRUE)
    c <- rnorm(N[m]*N[m],0,1)
    Probabilities <- matrix(exp(c/lambda)/(1+exp(c/lambda)),N[m])
    diag(Probabilities) <- 0
    prior[[m]]<-Probabilities
  }
  
  
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
  
  Kv <- 2*ncol(X[[1]]) + 1
  
  
  hyperparms <- list("network" = prior,
                     "theta0" = rep(0,Kv),
                     "invsigmatheta" = diag(Kv)/100,
                     "zeta0" = 0,
                     "insigma2zeta" = 2,
                     "a" = 4.2,
                     "b" = (4.2 - 2))
  
  # Set non informative prior for theta and zeta
  
  out <- peerMCMC(y, X, c(beta,gamma,0.4,se), hyperparms, iteration = MCMC.iteration)
  
  
  bayesian.estimate <- rbind(bayesian.estimate,apply(out[500:MCMC.iteration,], 2, mean))
  
  #save(bayesian.estimate, file = "bayesian.estimate.RData")
  
  out[[l]]          <- bayesian.estimate
}
t2 <- Sys.time()


t2 - t1
###########################
rm(list = ls())
setwd("/mnt/data/Thesis/ChapterI")
my.sum <- function(x) {
  out <- c(mean(x),
           sd(x),
           quantile(x, 0.25),
           median(x),
           quantile(x, 0.75))
  names(out) <- c("Mean", "Sd.", "1st Qu.", "Median", "3rd Qu.")
  return(out)
}

load("bayesian.estimate.RData")
(sum.bayesian.estimate <- t(apply(bayesian.estimate[1:200,],2,my.sum)))
write.csv(sum.bayesian.estimate, file = "sum.bayes.csv")

