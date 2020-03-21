#####################################Headline#####################################
rm(list = ls())
setwd("~/Dropbox/ARD/Simulations/CppFiles")
##################################################################################

####################################Simulations###################################
set.seed(123)

M <- 100 #Number of groups
N <- rep(50, M)

lambda <- c(0.5, 1, 1.5, 2)

Iteration <- 1000

save.cor <- matrix(NA, Iteration, length(lambda))

for (j in 1:length(lambda)) {
  for (i in 1:Iteration) {
    print(paste0("Iteration Lambda: ", j, "/", length(lambda), " - ", "Itetation MC: ", i, "/", Iteration))
    G <- list()
    X <- list()
    y <- list()
    
    
    ## True network
    prior <- list()
    
    G0 <- list()
    
    for (m in 1:M) {
      c <- rnorm(N[m] * N[m], 0, 1)
      Probabilities <- matrix(exp(c / lambda[j]) / (1 + exp(c / lambda[j])), N[m])
      diag(Probabilities) <- 0
      prior[[m]] <- Probabilities
      
      Gm <- matrix(runif(N[m] ^ 2), N[m], N[m]) < prior[[m]]
      diag(Gm) <- 0
      G0[[m]] <- Gm
    }
    
    
    
    # Estimated network
    G <- list()
    
    for (m in 1:M) {
      Gm <- matrix(runif(N[m] ^ 2), N[m], N[m]) < prior[[m]]
      diag(Gm) <- 0
      G[[m]] <- Gm
    }
    
    ## Cor G0 and G
    
    G0vec <- Gvec <- list()
    for (m in 1:M) {
      k <- 1:N[m]
      temp <- (k - 1) * N[m] + k
      G0vecm <- c(G0[[m]])
      G0vecm <- G0vecm[-temp]
      G0vec[[m]] <- G0vecm
      
      Gvecm <- c(G[[m]])
      Gvecm <- Gvecm[-temp]
      Gvec[[m]] <- Gvecm
    }
    
    G0vecall <- do.call("c", lapply(1:M, function(x) {
      G0vec[[x]]
    }))
    Gvecall <- do.call("c", lapply(1:M, function(x) {
      Gvec[[x]]
    }))
    
    
    save.cor[i, j] <- cor(G0vecall, Gvecall)
  }
}

my.sum <- function(x) {
  out <- c(mean(x),
           sd(x),
           quantile(x, 0.25),
           median(x),
           quantile(x, 0.75))
  names(out) <- c("Mean", "Sd.", "1st Qu.", "Median", "3rd Qu.")
  return(out)
}

apply(save.cor, 2, my.sum)

