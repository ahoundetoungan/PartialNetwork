#####################################Headline#####################################
rm(list = ls())
setwd("~/Dropbox/ARD/Simulations/Monte Carlo")
library(Rcpp)                         
library(AER)                       
sourceCpp("../CppFiles/Instruments.cpp")
sourceCpp("../CppFiles/DataForARD.cpp")
##################################################################################

####################################Simulations###################################
set.seed(123)

M <- 100 #Number of groups
N <- rep(50, M)

beta <- c(2, 1, 1.5)
gamma <- c(5,-3)
alpha <- 0.4
se <- 1

Iteration <- 1000

P <- matrix(1, N[1], N[1])
diag(P) <- 0
rs <- rowSums(P)
P <- P / rs
J <- diag(N[1]) - P


lest3    <- list()
lest3.1  <- matrix(NA, Iteration, 9)
lest3.2  <- matrix(NA, Iteration, 11)
colnames(lest3.1) <-
  c("Intercept",
    paste0("X", 1:2),
    paste0("GX", 1:2),
    "alpha",
    "Weak",
    "Wu",
    "Sargan")
colnames(lest3.2) <- 
  c(
    "Intercept",
    paste0("X", 1:2),
    paste0("GX", 1:2),
    paste0("GXc", 1:2),
    "alpha",
    "Weak",
    "Wu",
    "Sargan"
  )

for (l in 1:Iteration) {
  # Some useful variables
  W <- Y3 <- X <- GY3 <- GX <- GY3c <- GXc <- G2Xc <- list(M)
  
  #loop for group
  for (m in 1:M) {
    #Covariates
    X[[m]] <- cbind(rnorm(N[m], 0, 5), rpois(N[m], 6))
    
    feffect <-
      0.3 * X[[m]][1, 1] + 0.3 * X[[m]][3, 2] - 1.8 * X[[m]][50, 2]
    
    # Network depending on X
    Probabilities          <- matrix(0, N[m], N[m])
    for (i in 1:(N[m] - 1)) {
      for (j in (i + 1):N[m]) {
        Probabilities[i,j] <- pnorm(-4.5 + abs(X[[m]][i,1] - X[[m]][j,1]) - 2*abs(X[[m]][i,2] - X[[m]][j,2]))
      }
    }
    Probabilities          <- Probabilities + t(Probabilities)
    
    G <- Graph(Probabilities)
    #probtograph(Probabilities,G,gend,N)
    #True network row normalized
    rs <- rowSums(G)
    rs[rs == 1] <- 1
    W <- G / rs
    mean(rs)
    #True GX
    GX[[m]] <- W %*% X[[m]]
    
    #Y for section 3.3 Sub-Population Unobserved Fixed Effect
    Y3[[m]] <-
      solve(diag(rep(1, N[m])) - alpha * W) %*% (feffect + X[[m]] %*% beta[-1] +  GX[[m]] %*%
                                                   gamma + rnorm(N[m], 0, se))
    GY3[[m]] <- W %*% Y3[[m]]
    
    
    #Compute instruments
    # The function has 5 arguments
    # The network distribution distr and the explanatory variables X are required
    # y, S and pow are optional
    # S (by default = 2) is the number of draws of G needed for GX such that the first draw corresponds to the draw G used to compute GY
    # pow is used to compute GX, G^2X, ..., G^(pow)X
    
    # The output of this function is as follow
    # $GY  estimates GY by drawin one G from the network distribution and multiplies it by Y
    # $GX is a list of S instrument matrices, Each index in the list corresponds to one draw of G to compute GX. The first draw is also used to
    # compute GY
    # Each matrix instrument is cube (three dimensions) (N,K,pow)
    # the first axe indexes the individual, the second axe index the variable,
    # and the fird axe the power (GX, G^2X, G^3X)
    # For example $GX[[2]][,,3] is G^3X where G is the second draw of G
    
    instr3 <- instruments(Probabilities, X[[m]], Y3[[m]], pow = 2)
    
    GY3c[[m]] <- instr3$GY
    
    GXc[[m]] <- instr3$GX[[1]][, , 1]
    G2Xc[[m]] <- instr3$GX[[2]][, , 2]
    #print(paste("Iteration:",l,"** Group:",m))
  }
  
  
  # Concatenate M groups data
  
  # Y for section 3.1, 3.2 and 3.3
  Y3allm0 <- do.call("c", lapply(1:M, function(x) {
    J %*% Y3[[x]]
  }))
  
  # In section 3.1 or 3.3 GY may be observed
  GY3all <- do.call("c", lapply(1:M, function(x) {
    J %*% GY3[[x]]
  }))
  
  # GY constructed for section 3.3
  GY3callm0 <- do.call("c", lapply(1:M, function(x) {
    J %*% GY3c[[x]]
  }))
  
  #Centered covariates
  Xallm0 <- do.call(rbind, lapply(1:M, function(x) {
    J %*% X[[x]]
  }))
  
  # GX
  GXallm0 <- do.call(rbind, lapply(1:M, function(x) {
    J %*% GX[[x]]
  }))
  
  # G^pX constructed
  GXcallm0 <-
    do.call(rbind, lapply(1:M, function(x) {
      J %*% GXc[[x]]
    }))
  G2Xcallm0 <-
    do.call(rbind, lapply(1:M, function(x) {
      J %*% G2Xc[[x]]
    }))
  
  ###### Estimation
  # estimation section 3.3
  est3.1 <-
    ivreg(Y3allm0 ~ Xallm0 + GXallm0 + GY3all |
            Xallm0 + GXallm0 + G2Xcallm0)
  sest3.1 <- summary(est3.1, diagnostic = TRUE)
  lest3.1[l,] <-
    c(sest3.1$coefficients[, 1], sest3.1$diagnostics[, 3])
  #print(paste("Iteration:", l, "Model fixed effects, Gy is observed"))
  #print(lest3.1[l,])
  
  #Y is not observed
  est3.2 <-
    ivreg(
      Y3allm0 ~ Xallm0 + GXallm0 + GXcallm0 + GY3callm0 |
        Xallm0 + GXallm0 + GXcallm0 + G2Xcallm0
    )
  sest3.2 <- summary(est3.2, diagnostic = TRUE)
  lest3.2[l,] <-
    c(sest3.2$coefficients[, 1], sest3.2$diagnostics[, 3])
  #print(paste0("lambda: ", i,"/",length(lambda),  " - Iteration: ", l, "/", Iteration))
  print(paste0("Iteration: ", l, "/", Iteration))
  #print(lest3.2[l,])
}

#Lest3[[i]] <- list(lest3.1, lest3.2)

save(lest3.1, lest3.2, file = "Results/Lest3x.Rdata")

################################################
load("Results/Lest3x.Rdata")


my.sum <- function(x) {
  out <- c(mean(x),
           sd(x),
           quantile(x, 0.25),
           median(x),
           quantile(x, 0.75))
  names(out) <- c("Mean", "Sd.", "1st Qu.", "Median", "3rd Qu.")
  return(out)
}


(sum.lest3.1 <- round(t(as.matrix(apply(lest3.1 , 2, my.sum))),3))
(sum.lest3.2 <- round(t(as.matrix(apply(lest3.2 , 2, my.sum))),3))

write.csv(sum.lest3.1, file = "Results/sum.lest3.1x.csv")
write.csv(sum.lest3.2, file = "Results/sum.lest3.2x.csv")
