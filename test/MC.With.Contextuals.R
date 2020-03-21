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
lambda <- c(seq(0.01, 2, 0.01), Inf)

beta <- c(2, 1, 1.5)
gamma <- c(5,-3)
alpha <- 0.4
se <- 1

Iteration <- 1000

Lest2 <- list()
for (i in 1:length(lambda)) {
  lest2.1 <- matrix(NA, Iteration, 9)
  lest2.2 <- matrix(NA, Iteration, 11)
  
  colnames(lest2.1) <-
    c(
      "Intercept",
      paste0("X", 1:2),
      paste0("GX", 1:2),
      "alpha",
      "Weak",
      "Wu",
      "Sargan"
    )
  colnames(lest2.2) <-
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
    W <- Y2 <- X <- GY2 <- GX <- GY2c <- GXc <- G2Xc <- list(M)
    
    #loop for group
    for (m in 1:M) {
      c <- rnorm(N[m] * N[m], 0, 1)
      
      Probabilities <-
        matrix(exp(c / lambda[i]) / (1 + exp(c / lambda[i])), N[m])
      diag(Probabilities) <- 0
      
      G <- Graph(Probabilities)
      #probtograph(Probabilities,G,gend,N)
      #True network row normalized
      rs <- rowSums(G)
      rs[rs == 1] <- 1
      W <- G / rs
      
      #Covariates
      X[[m]] <- cbind(rnorm(N[m], 0, 5), rpois(N[m], 6))
      
      #True GX
      GX[[m]] <- W %*% X[[m]]
      
      #Y for section 3.2 with contextual effects
      Y2[[m]] <-
        solve(diag(rep(1, N[m])) - alpha * W) %*% (cbind(rep(1, N[m]), X[[m]]) %*%
                                                     beta +  GX[[m]] %*% gamma + rnorm(N[m], 0, se))
      GY2[[m]] <- W %*% Y2[[m]]
      
      
      
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
      
      instr2 <- instruments(Probabilities, X[[m]], Y2[[m]], pow = 2)
      
      GY2c[[m]] <- instr2$GY
      
      GXc[[m]] <- instr2$GX[[1]][, , 1]
      G2Xc[[m]] <- instr2$GX[[2]][, , 2]
      #print(paste("Iteration:",l,"** Group:",m))
    }
    
    
    # Concatenate M groups data
    
    # Y for section 3.1, 3.2 and 3.3
    Y2all <- do.call("c", lapply(1:M, function(x) {
      Y2[[x]]
    }))
    
    # In section 3.1 or 3.3 GY may be observed
    GY2all <- do.call("c", lapply(1:M, function(x) {
      GY2[[x]]
    }))
    
    # GY constructed for section 3.2
    GY2call <- do.call("c", lapply(1:M, function(x) {
      GY2c[[x]]
    }))
    
    
    # Covariates
    Xall <- do.call(rbind, lapply(1:M, function(x) {
      X[[x]]
    }))
    
    # GX
    GXall <- do.call(rbind, lapply(1:M, function(x) {
      GX[[x]]
    }))
    
    # G^pX constructed
    GXcall <- do.call(rbind, lapply(1:M, function(x) {
      GXc[[x]]
    }))
    G2Xcall <- do.call(rbind, lapply(1:M, function(x) {
      G2Xc[[x]]
    }))
    
    ###### Estimation
    # estimation section 3.2
    #GY is observed
    est2.1 <-
      ivreg(Y2all ~ Xall + GXall + GY2all |
              Xall  + GXall + G2Xcall)
    sest2.1 <- summary(est2.1, diagnostic = TRUE)
    lest2.1[l, ] <-
      c(sest2.1$coefficients[, 1], sest2.1$diagnostics[, 3])
    #print(paste("Iteration:", l, "Model with contextual effects"))
    #print(lest2.1[l, ])
    
    #GY is not observed
    est2.2 <-
      ivreg(Y2all ~ Xall + GXall +  GXcall + GY2call |
              Xall  + GXall + GXcall + G2Xcall)
    sest2.2 <- summary(est2.2, diagnostic = TRUE)
    lest2.2[l, ] <-
      c(sest2.2$coefficients[, 1], sest2.2$diagnostics[, 3])
    print(paste0("lambda: ", i,"/",length(lambda),  " - Iteration: ", l, "/", Iteration))
    #print(lest2.2[l, ])
  }
  
  
  Lest2[[i]] <- list(lest2.1, lest2.2)
}


save(Lest2, lambda, file = "Results/Lest2.Rdata")

################################################
library(ggplot2)
multiplot <-
  function(...,
           plotlist = NULL,
           file,
           cols = 1,
           layout = NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                       ncol = cols,
                       nrow = ceiling(numPlots / cols))
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
        
        print(plots[[i]],
              vp = viewport(
                layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col
              ))
      }
    }
  }

################################################
rm(list = ls())
load("Results/Lest2.Rdata")

my.sum <- function(x) {
  out <- c(mean(x),
           sd(x),
           quantile(x, 0.25),
           median(x),
           quantile(x, 0.75))
  names(out) <- c("Mean", "Sd.", "1st Qu.", "Median", "3rd Qu.")
  return(out)
}



lindex  <- which(lambda %in%  seq(0.01, 2, 0.01))
lambdai <- Inf
temp <- Lest2[[which(round(lambda, 2) == lambdai)]]

lambda[which(round(lambda, 2) == lambdai)]

lest2.1 <- temp[[1]]
lest2.2 <- temp[[2]]

(sum.lest2.1 <- round(t(as.matrix(apply(lest2.1 , 2, my.sum))),3))
(sum.lest2.2 <- round(t(as.matrix(apply(lest2.2 , 2, my.sum))),3))

write.csv(sum.lest2.1, file = paste("Results/sum.lest2.1_", lambdai, ".csv"))
write.csv(sum.lest2.2, file = paste("Results/sum.lest2.2_", lambdai, ".csv"))


#### Graph
alpha.mean.2.1 <-
  unlist(lapply(lindex, function(w) {
    mean(Lest2[[w]][[1]][, "alpha"])
  }))
alpha.1st.Qu.2.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest2[[w]][[1]][, "alpha"], 0.25)
  }))
alpha.2st.Qu.2.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest2[[w]][[1]][, "alpha"], 0.75)
  }))
alpha.INF.CI.2.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest2[[w]][[1]][, "alpha"], 0.025)
  }))
alpha.SUP.CI.2.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest2[[w]][[1]][, "alpha"], 0.975)
  }))

tab.graph.2.1 <-
  data.frame(
    lambda = lambda[lindex],
    mean = alpha.mean.2.1,
    Inf.IQI = alpha.1st.Qu.2.1,
    Sup.IQI = alpha.2st.Qu.2.1,
    Inf.CI = alpha.INF.CI.2.1,
    Sup.CI = alpha.SUP.CI.2.1
  )

alpha.mean.2.2 <-
  unlist(lapply(lindex, function(w) {
    mean(Lest2[[w]][[2]][, "alpha"])
  }))
alpha.1st.Qu.2.2 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest2[[w]][[2]][, "alpha"], 0.25)
  }))
alpha.2st.Qu.2.2 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest2[[w]][[2]][, "alpha"], 0.75)
  }))
alpha.INF.CI.2.2 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest2[[w]][[2]][, "alpha"], 0.025)
  }))
alpha.SUP.CI.2.2 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest2[[w]][[2]][, "alpha"], 0.975)
  }))

tab.graph.2.2 <-
  data.frame(
    lambda = lambda[lindex],
    mean = alpha.mean.2.2,
    Inf.IQI = alpha.1st.Qu.2.2,
    Sup.IQI = alpha.2st.Qu.2.2,
    Inf.CI = alpha.INF.CI.2.2,
    Sup.CI = alpha.SUP.CI.2.2
  )



g1 <- ggplot(tab.graph.2.1, aes(x = lambda)) + theme_bw() +
  geom_ribbon(aes(ymin = Inf.CI, ymax = Sup.CI, fill = "95% CI"), alpha = 1) +
  geom_ribbon(aes(ymin = Inf.IQI, ymax = Sup.IQI, fill = "Interquartile Interval"), alpha = 1) +
  geom_line(aes(y = mean, colour = "Mean")) + ylab("") +
  scale_fill_manual(NULL, values = c("#ccccff", "#55ffff")) +
  scale_colour_manual(NULL, values = "#ff0000") +
  theme(
    legend.position = c(0.3, 0.88),
    legend.box = 'horizontal',
    axis.ticks = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 1)) +
  ylab("alpha") +
  ggtitle(paste(expression(Gy), "is oberved")) +
  scale_y_continuous(breaks = seq(0.36, 0.425, 0.02),
                     limits = c(0.375, 0.425))+
  scale_x_continuous(breaks = seq(0,2,0.5), labels = as.expression(seq(0,2,0.5)), 
                     sec.axis = sec_axis(~ ., name = "Correlation",
                                         labels = c(1,0.394,0.173,0.092,0.056)))

g2 <-
  ggplot(tab.graph.2.2, aes(x = lambda)) +  theme_bw() +
  geom_ribbon(aes(ymin = Inf.CI, ymax = Sup.CI, fill = "95% CI"), alpha = 1) +
  geom_ribbon(aes(ymin = Inf.IQI, ymax = Sup.IQI, fill = "Interquartile Interval"), alpha = 1) +
  geom_line(aes(y = mean, colour = "Mean")) + ylab("") +
  scale_fill_manual(NULL, values = c("#ccccff", "#55ffff")) +
  scale_colour_manual(NULL, values = "#ff0000") +
  theme(
    legend.position = c(0.3, 0.88),
    legend.box = 'horizontal',
    axis.ticks = element_blank()
  ) +
  ylab("") +
  ggtitle(paste(expression(Gy), "is not oberved")) +
  scale_y_continuous(breaks = seq(0.36, 0.425, 0.02),
                     limits = c(0.375, 0.425)) +
  scale_x_continuous(breaks = seq(0,2,0.5), labels = as.expression(seq(0,2,0.5)), 
                     sec.axis = sec_axis(~ ., name = "Correlation",
                                         labels = c(1,0.394,0.173,0.092,0.056)))

multiplot(g1, g2, cols = 2)
# size 11x4.5 inch