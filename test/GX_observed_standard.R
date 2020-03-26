##################################################################################
# This code replicates the Monte Carlo simulations when GX is observed and       #
# the precision of simulated networks are controled by lambda.                   #
# We distinguish the case where Gy is observed and not observed (section 3)      #
#####################################Headline#####################################
rm(list = ls())
library(AER)                          # Install AER if not already done.
library(PartialNetwork)               # Install PartialNetwork if not already done.
library(doParallel)                   # To run the Monte Carlo in parallel
##################################################################################
# our summary function
our.sum <- function(x) {
  out <- c(mean(x, na.rm = TRUE),
           sd(x, na.rm = TRUE),
           quantile(x, c(0.025, 0.25), na.rm = TRUE),
           median(x, na.rm = TRUE),
           quantile(x, c(0.75, 0.975), na.rm = TRUE))
  names(out) <- c("Mean", "Sd.", "pct2.5", "1st Qu.", "Median", "3rd Qu.", "pct97.5")
  return(out)
}


# function to perform the simulation
# l stands for the l-th simulation
# lambda network precision parameter
fsim <- function(l, lambda){
  cat("lambda ", lambda, " -- Iteration ", l, "\n")
  M          <- 100          # Number of groups
  N          <- rep(50,M)   # Group size
  # Parameters
  
  beta      <- c(2,1,1.5)
  gamma     <- c(5,-3)
  alpha     <- 0.4
  se        <- 1
  
  P1        <- matrix(1, N[1], N[1])
  diag(P1)  <- 0
  rs        <- rowSums(P1)
  P1        <- P1 / rs
  J         <- diag(N[1]) - P1
  
  
  # Some useful variables
  W     <- Y1    <- Y2   <- Y3    <- X     <- GY1   <- GY2    <- GY3      <- GX  <- list(M)
  GY1c  <- GY2c  <- GY3c <- GX1c  <- G2X1c <-GX1c0  <- G2X1c0 <- indexall <- list(M)
  GX2c  <- G2X2c <-GX2c0 <-G2X2c0 <- GX3c  <- G2X3c <- GX3c0  <-G2X3c0    <- list(M)
  
  Gvec  <- c() 
  G0vec <- c()
  
  
  #loop over groups
  for (m in 1:M) {
    #Generate link probabilities
    c                <- rnorm(N[m] * N[m], 0, 1)
    
    distr            <- matrix(1 / (1 + exp(-c / lambda)), N[m])
    diag(distr)      <- 0
    
    #The complete graph
    G                <- sim.network(dnetwork = distr)
    
    #To mesure the network precision we sample another graph
    G0               <- sim.network(dnetwork = distr)
    #store off diagonal elements of G and G0
    #we will compute the graph correlation
    Gvec             <- c(Gvec, G[row(G) != col(G)])
    G0vec            <- c(G0vec, G0[row(G0) != col(G0)])
    
    #True network row normalized
    rs               <- rowSums(G)
    rs[rs == 0]      <- 1
    W[[m]]           <- G / rs
    
    #Covariates
    X[[m]]   <- cbind(rnorm(N[m],0,5),rpois(N[m],6))
    
    #Fixed effects
    feffect  <- 0.3 * X[[m]][1, 1] + 0.3 * X[[m]][3, 2] - 1.8 * X[[m]][50, 2]
    
    #True GX
    GX[[m]]  <- W[[m]] %*% X[[m]]
    
    #Y for the section without contextual effects
    Y1[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1,N[m]), X[[m]])%*%beta + rnorm(N[m],0,se))
    GY1[[m]] <- W[[m]] %*% Y1[[m]]
    
    #Y for thesection with contextual effects
    Y2[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (cbind(rep(1, N[m]), X[[m]]) %*% beta +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY2[[m]] <- W[[m]] %*% Y2[[m]]
    
    #Y for the section with Sub-Population Unobserved Fixed Effect
    Y3[[m]]  <- solve(diag(rep(1, N[m])) - alpha * W[[m]]) %*% (feffect + X[[m]] %*% beta[-1] +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
    GY3[[m]] <- W[[m]] %*% Y3[[m]]
    
    #Compute instruments
    instr1 <- sim.IV(dnetwork = distr, X[[m]], Y1[[m]], replication = 1, power = 2)
    instr2 <- sim.IV(dnetwork = distr, X[[m]], Y2[[m]], replication = 1, power = 2)
    instr3 <- sim.IV(dnetwork = distr, X[[m]], Y3[[m]], replication = 1, power = 2)
    
    GY1c[[m]]   <- instr1[[1]]$G1y
    GX1c0[[m]]  <- instr1[[1]]$G1X[, , 1] 
    G2X1c0[[m]] <- instr1[[1]]$G1X[, , 2]  
    GX1c[[m]]   <- instr1[[1]]$G2X[, , 1]
    G2X1c[[m]]  <- instr1[[1]]$G2X[, , 2]
    
    
    GY2c[[m]]   <- instr2[[1]]$G1y
    GX2c0[[m]]  <- instr2[[1]]$G1X[, , 1] 
    G2X2c0[[m]] <- instr2[[1]]$G1X[, , 2]  
    GX2c[[m]]   <- instr2[[1]]$G2X[, , 1]
    G2X2c[[m]]  <- instr2[[1]]$G2X[, , 2]
    
    
    GY3c[[m]]   <- instr3[[1]]$G1y
    GX3c0[[m]]  <- instr3[[1]]$G1X[, , 1] 
    G2X3c0[[m]] <- instr3[[1]]$G1X[, , 2]  
    GX3c[[m]]   <- instr3[[1]]$G2X[, , 1]
    G2X3c[[m]]  <- instr3[[1]]$G2X[, , 2]
  }
  
  # Concatenate M groups data
  # Y 
  Y1all     <- do.call("c", lapply(1:M, function(x) Y1[[x]]))
  Y2all     <- do.call("c", lapply(1:M, function(x) Y2[[x]]))
  Y3all     <- do.call("c", lapply(1:M, function(x) J %*% Y3[[x]]))
  
  
  # GY observed
  GY1all    <- do.call("c", lapply(1:M, function(x) GY1[[x]]))
  GY2all    <- do.call("c", lapply(1:M, function(x) GY2[[x]]))
  GY3all    <- do.call("c", lapply(1:M, function(x) J %*% GY3[[x]]))
  
  # GY constructed 
  GY1call   <- do.call("c", lapply(1:M, function(x) GY1c[[x]]))
  GY2call   <- do.call("c", lapply(1:M, function(x) GY2c[[x]]))
  GY3call   <- do.call("c", lapply(1:M, function(x) J %*% GY3c[[x]]))
  
  
  # X
  Xall      <- do.call(rbind, lapply(1:M, function(x) X[[x]]))
  X3all     <- do.call(rbind, lapply(1:M, function(x) J %*% X[[x]]))
  
  # GX observed
  GXall     <- do.call(rbind, lapply(1:M, function(x) GX[[x]]))
  GX3all    <- do.call(rbind, lapply(1:M, function(x) J %*% GX[[x]]))
  
  
  # G^pX constructed
  # model without contextual effects
  # same draw
  GX1c0all   <- do.call(rbind, lapply(1:M, function(x) GX1c0[[x]]))
  G2X1c0all  <- do.call(rbind, lapply(1:M, function(x) G2X1c0[[x]]))
  # different draw
  GX1call    <- do.call(rbind, lapply(1:M, function(x) GX1c[[x]]))
  G2X1call   <- do.call(rbind, lapply(1:M, function(x) G2X1c[[x]]))
  
  # model with contextual effects
  # same draw
  GX2c0all   <- do.call(rbind, lapply(1:M, function(x) GX2c0[[x]]))
  G2X2c0all  <- do.call(rbind, lapply(1:M, function(x) G2X2c0[[x]]))
  # different draw
  GX2call    <- do.call(rbind, lapply(1:M, function(x) GX2c[[x]]))
  G2X2call   <- do.call(rbind, lapply(1:M, function(x) G2X2c[[x]]))

  # model with fixed effects
  # same draw
  GX3c0all   <- do.call(rbind, lapply(1:M, function(x) J %*% GX3c0[[x]]))
  G2X3c0all  <- do.call(rbind, lapply(1:M, function(x) J %*% G2X3c0[[x]]))
  # different draw
  GX3call    <- do.call(rbind, lapply(1:M, function(x) J %*% GX3c[[x]]))
  G2X3call   <- do.call(rbind, lapply(1:M, function(x) J %*% G2X3c[[x]]))
  

  ###### Estimation
  #estimation section 3.1
  #if Gy is observed. We use GX constructed as instruments
  sest1.1.1   <- summary(ivreg(Y1all ~ Xall + GY1all | Xall +  GX1call), diagnostic=TRUE)
  lest1.1.1   <- c(sest1.1.1$coefficients[,1],sest1.1.1$diagnostics[,3])
  
  #if Gy is observed. We use GX and GGX constructed as instruments
  sest1.1.2   <- summary(ivreg(Y1all ~ Xall + GY1all | Xall +  GX1call + G2X1call) , diagnostic=TRUE)
  lest1.1.2   <- c(sest1.1.2$coefficients[,1],sest1.1.2$diagnostics[,3])
  
  #if Gy is not observed. 
  #Same draw
  #We use GX constructed as instruments
  sest1.2.1.1 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1c0all) , diagnostic=TRUE)
  lest1.2.1.1 <- c(sest1.2.1.1$coefficients[,1],sest1.2.1.1$diagnostics[,3], cor((GY1all-GY1call),GX1c0all))
  
  #We use GX and GGX constructed as instruments
  sest1.2.1.2 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1c0all + G2X1c0all), diagnostic=TRUE)
  lest1.2.1.2 <- c(sest1.2.1.2$coefficients[,1],sest1.2.1.2$diagnostics[,3], cor((GY1all-GY1call),cbind(GX1c0all,G2X1call)))
  
  #different draw
  #We use GX constructed as instruments
  sest1.2.2.1 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1call), diagnostic=TRUE)
  lest1.2.2.1 <- c(sest1.2.2.1$coefficients[,1],sest1.2.2.1$diagnostics[,3], cor((GY1all-GY1call),GX1call))
  
  #We use GX and GGX constructed as instruments
  sest1.2.2.2 <- summary(ivreg(Y1all ~ Xall + GY1call | Xall +  GX1call + G2X1call), diagnostic=TRUE)
  lest1.2.2.2 <- c(sest1.2.2.2$coefficients[,1],sest1.2.2.2$diagnostics[,3], cor((GY1all-GY1call),cbind(GX1call,G2X1call)))
  
  #estimation section 3.2
  #GY is observed
  sest2.1     <- summary(ivreg(Y2all ~ Xall + GXall + GY2all | Xall  + GXall + G2X2call), diagnostic = TRUE)
  lest2.1     <- c(sest2.1$coefficients[, 1], sest2.1$diagnostics[, 3])
  
  #GY is not observed
  sest2.2     <- summary(ivreg(Y2all ~ Xall + GXall +  GX2c0all + GY2call | Xall  + GXall + GX2c0all + G2X2call), diagnostic = TRUE)
  lest2.2     <- c(sest2.2$coefficients[, 1], sest2.2$diagnostics[, 3])
  if (lambda == 0) {
    lest2.2   <- c(sest2.2$coefficients[1:5, 1], c(NA, NA), sest2.2$coefficients[6, 1], sest2.2$diagnostics[, 3])
  }
  
  # estimation section 3.3
  #Gy is observed
  sest3.1     <- summary(ivreg(Y3all ~ X3all  + GX3all + GY3all | X3all + GX3all + G2X3call), diagnostic = TRUE)
  lest3.1     <- c(sest3.1$coefficients[, 1], sest3.1$diagnostics[, 3])
  
  #Gy is not observed
  sest3.2     <- summary(ivreg(Y3all ~ X3all + GX3all + GX3c0all + GY3call | X3all + GX3all + GX3c0all + G2X3call), diagnostic = TRUE)
  lest3.2     <- c(sest3.2$coefficients[, 1], sest3.2$diagnostics[, 3])
  if (lambda == 0) {
    lest3.2   <- c(sest3.2$coefficients[1:5, 1], c(NA, NA), sest3.2$coefficients[6, 1], sest3.2$diagnostics[, 3])
  }
  
  #network correlation
  cornet      <- cor(Gvec, G0vec)
  
  c(lest1.1.1, lest1.1.2, lest1.2.1.1, lest1.2.1.2, lest1.2.2.1,
    lest1.2.2.2, lest2.1, lest2.2, lest3.1, lest3.2, cornet)
}

# monte carlo function for each parlambda
f.mc <- function(iteration, lambda) {
  out.mc        <- mclapply(1:iteration, function(w) fsim(w, lambda), mc.cores = 8L)

  simu          <- t(do.call(cbind, out.mc))
  
  # the colnames
  tmp <- c("Intercept",paste0("X",1:2),"alpha","Weak","Wu","Sargan")
  c1  <- paste0("No Con - GY obs - ins GX ", tmp)
  c2  <- paste0("No Con - GY obs - ins GX GGX ", tmp)
  c3  <- paste0("No Con - GY notobs - ins GX - sam draw ", c(tmp,"corGX1e","corGX2e"))
  c4  <- paste0("No Con - GY notobs - ins GX GGX - sam draw ", c(tmp,"corGX1e","corGX2e","corGGX1e","corGGX2e"))
  c5  <- paste0("No Con - GY notobs - ins GX - dif draw ", c(tmp,"corGX1e","corGX2e"))
  c6  <- paste0("No Con - GY notobs - ins GX GGX - dif draw ", c(tmp,"corGX1e","corGX2e","corGGX1e","corGGX2e"))
  
  tmp <- c("Intercept", paste0("X", 1:2), paste0("GX", 1:2), "alpha", "Weak", "Wu", "Sargan")
  c7  <- paste0("Wit Con - GY obs ", tmp)
  c9  <- paste0("Fix eff - GY obs ", tmp)
  tmp <- c("Intercept", paste0("X", 1:2), paste0("GX", 1:2), paste0("GXc", 1:2), "alpha", "Weak", "Wu", "Sargan")
  c8  <- paste0("Wit Con - GY notobs ", tmp)
  c10 <- paste0("Fix eff - GY notobs ", tmp)
  c11 <- "Network correlation"
  
  colnames(simu) <- c(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11)
  
  # summary for all simulation using ARD
  results        <- t(apply(simu, 2, our.sum))
  results
}

set.seed(123)

# Number of simulations
iteration <- 1000
# Monte Carlo for several values of lambda
veclambda <- c(seq(0, 2, 0.01), Inf)
out       <- lapply(veclambda, function(lambda) f.mc(iteration, lambda))

# results for specific lambda
lambda    <- 1
ilambda   <- which(veclambda == lambda)
out[[ilambda]]

lambda    <- Inf
ilambda   <- which(veclambda == lambda)
out[[ilambda]]

############## plot estimation of alpha for lambda from 0 to 2 ##################

### multiplot function
library(ggplot2)
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots    <- c(list(...), plotlist)
  
  numPlots <- length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)), ncol = cols, nrow = ceiling(numPlots / cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}
row.n     <- rownames(out[[1]])

# fonction to create graph data
# depends on row name of alpha in out
f.data    <- function(row.n.i) {
  tmp           <- lapply(out[-202], function(w) c(w[row.n == row.n.i, c(1, 3, 7, 4, 6)], w[95,1]))
  tmp           <- t(do.call(cbind, tmp))
  colnames(tmp) <- c("Mean", "pct2.5", "pct97.5", "pct25", "pct75", "correlation")
  as.data.frame(cbind("lambda" = veclambda[-202], tmp))
}

# function for the graph 
f.graph    <- function(row.n.i) {
  
  dat      <- f.data(row.n.i)
  breaks   <- seq(0.35, 0.45, 0.02)
  limits   <- c(0.35, 0.45)
  if (grepl("No Con", row.n.i)) {
    breaks <- seq(0.36, 0.425, 0.02)
    limits <- c(0.375, 0.425)
  }
  if (grepl("Fix eff", row.n.i)) {
    breaks <- seq(-6.8, 6.7, 1.8)
    limits <- c(-7.5, 7.5)
  }
    ggplot(dat, aes(x = lambda)) + theme_bw()  +
    geom_ribbon(aes(ymin = pct2.5, ymax = pct97.5, fill = "95% CI"), alpha = 1) +
    geom_ribbon(aes(ymin = pct25, ymax = pct75, fill = "Interquartile Interval"), alpha = 1) +
    geom_line(aes(y = Mean, colour = "Mean")) + ylab("") + scale_fill_manual(NULL, values = c("#ccccff", "#55ffff")) +
    scale_colour_manual(NULL, values = "#ff0000") +
    theme(legend.position = c(0.3, 0.885), legend.box = 'horizontal', axis.ticks = element_blank()) +
    guides(colour = guide_legend(nrow = 1)) + ylab("alpha") +
    ggtitle(paste(expression(Gy), ifelse(grepl("GY obs", row.n.i), "is oberved", "is not oberved"))) +
    scale_y_continuous(breaks = breaks, limits = limits) +
    scale_x_continuous(breaks = seq(0,2,0.5), labels = as.expression(seq(0,2,0.5)), 
                       sec.axis = sec_axis(~ ., name = "Correlation", labels = round(dat$correlation[dat$lambda %in% seq(0,2,0.5)], 3)))
}

# model without contextual effects
g11 <- f.graph("No Con - GY obs - ins GX alpha")
g12 <- f.graph("No Con - GY notobs - ins GX - dif draw alpha")
multiplot(g11, g12, cols = 2)
# size 10x4.35 inch

# model with contextual effects
g21 <- f.graph("Wit Con - GY obs alpha")
g22 <- f.graph("Wit Con - GY notobs alpha")
multiplot(g21, g22, cols = 2)
# size 10x4.35 inch

# model with fixed effects
g31 <- f.graph("Fix eff - GY obs alpha" )
g32 <- f.graph("Fix eff - GY notobs alpha" )
multiplot(g31, g32, cols = 2)
# size 10x4.35 inch