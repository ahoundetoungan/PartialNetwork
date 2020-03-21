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

Lest1 <- list()
for (i in 1:length(lambda)) {
  lest1.1.1 <- lest1.1.2 <- matrix(NA,Iteration,7)
  lest1.2.1.1 <- lest1.2.2.1 <- matrix(NA,Iteration,9)
  lest1.2.1.2 <- lest1.2.2.2 <- matrix(NA,Iteration,11)
  
  colnames(lest1.1.1) <- colnames(lest1.1.2) <- c("Intercept",paste0("X",1:2),"alpha","Weak","Wu","Sargan")
  colnames(lest1.2.1.1) <- colnames(lest1.2.2.1) <- c(colnames(lest1.1.2),"corGX1e","corGX2e")
  colnames(lest1.2.1.2) <- colnames(lest1.2.2.2) <- c(colnames(lest1.1.2),"corGX1e","corGX2e","corGGX1e","corGGX2e")
  
  for (l in 1:Iteration) {
    #print(paste0("iteration: ", l, "/", Iteration))
    # Some useful variables
    W<-Y1<-X<-GY1<-GX<-GY1c<-GXc0<-G2Xc0<-GXc<-G2Xc<-list(M)
    
    #loop for group
    for (m in 1:M) {
      
      c <- rnorm(N[m]*N[m],0,1)
      
      Probabilities <- matrix(exp(c / lambda[i]) / (1 + exp(c / lambda[i])), N[m])
      diag(Probabilities) <- 0
      
      G<-Graph(Probabilities)
      #probtograph(Probabilities,G,gend,N)
      #True network row normalized
      rs<-rowSums(G)
      rs[rs==1]<-1
      W<-G/rs
      
      #Covariates
      X[[m]]<-cbind(rnorm(N[m],0,5),rpois(N[m],6))
      
      #True GX
      GX[[m]]<-W%*%X[[m]]
      
      #Y for sections 3.1 without contextual effects
      Y1[[m]]<-solve(diag(rep(1,N[m]))-alpha*W)%*%(cbind(rep(1,N[m]), X[[m]])%*%beta + rnorm(N[m],0,se))
      GY1[[m]]<-W%*%Y1[[m]]
      
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
      
      instr1<-instruments(Probabilities, X[[m]], Y1[[m]], pow=2)
      
      GY1c[[m]] <- instr1$GY
      
      GXc0[[m]]<-instr1$GX[[1]][,,1] 
      G2Xc0[[m]]<-instr1$GX[[1]][,,2] 
      GXc[[m]]<-instr1$GX[[2]][,,1] 
      G2Xc[[m]]<-instr1$GX[[2]][,,2]
      #print(paste("Iteration:",l,"** Group:",m))
    }
    
    
    # Concatenate M groups data 
    
    # Y for section 3.1, 3.2 and 3.3 
    Y1all<-do.call("c",lapply(1:M, function(x){Y1[[x]]}))
    
    # In section 3.1 or 3.3 GY may be observed
    GY1all<-do.call("c",lapply(1:M, function(x){GY1[[x]]}))
    
    # GY constructed for section 3.1
    GY1call<-do.call("c",lapply(1:M, function(x){GY1c[[x]]}))
    
    
    # Covariates
    Xall<-do.call(rbind,lapply(1:M, function(x){X[[x]]}))
    
    # GX
    GXall<-do.call(rbind,lapply(1:M, function(x){GX[[x]]}))
    
    # G^pX constructed
    GXc0all<-do.call(rbind,lapply(1:M, function(x){GXc0[[x]]}))
    GXcall<-do.call(rbind,lapply(1:M, function(x){GXc[[x]]}))
    G2Xcall<-do.call(rbind,lapply(1:M, function(x){G2Xc[[x]]}))
    G2Xc0all<-do.call(rbind,lapply(1:M, function(x){G2Xc0[[x]]}))
    
    ###### Estimation
    # estimation section 3.1
    # if Gy is observed. We use GX constructed as instruments
    est1.1.1<-ivreg(Y1all ~ Xall + GY1all | Xall +  GXcall)
    sest1.1.1<-summary(est1.1.1, diagnostic=TRUE)
    lest1.1.1[l,] <-  c(sest1.1.1$coefficients[,1],sest1.1.1$diagnostics[,3])
    #print(paste("Iteration:",l,"Model without contextual effects, Gy observed and GX is instrument"))
    #print(lest1.1.1[l,])
    
    # if Gy is observed. We use GX and GGX constructed as instruments
    est1.1.2<-ivreg(Y1all ~ Xall + GY1all | Xall +  GXcall + G2Xcall)
    sest1.1.2<-summary(est1.1.2, diagnostic=TRUE)
    lest1.1.2[l,] <-  c(sest1.1.2$coefficients[,1],sest1.1.2$diagnostics[,3])
    #print(paste("Iteration:",l,"Model without contextual effects, Gy observed GX, GGX are instruments"))
    #print(lest1.1.2[l,])
    
    # if Gy is not observed. 
    #Same draw
    #We use GX constructed as instruments
    est1.2.1.1<-ivreg(Y1all ~ Xall + GY1call | Xall +  GXc0all)
    sest1.2.1.1<-summary(est1.2.1.1, diagnostic=TRUE)
    lest1.2.1.1[l,] <- c(sest1.2.1.1$coefficients[,1],sest1.2.1.1$diagnostics[,3],cor((GY1all-GY1call),GXc0all))
    #print(paste("Iteration:",l,"Model without contextual effects, Gy not observed, same draw and GX is instrument"))
    #print(lest1.2.1.1[l,])
    
    #We use GX and GGX constructed as instruments
    est1.2.1.2<-ivreg(Y1all ~ Xall + GY1call | Xall +  GXc0all + G2Xc0all)
    sest1.2.1.2<-summary(est1.2.1.2, diagnostic=TRUE)
    lest1.2.1.2[l,] <- c(sest1.2.1.2$coefficients[,1],sest1.2.1.2$diagnostics[,3],cor((GY1all-GY1call),cbind(GXc0all,G2Xcall)))
    #print(paste("Iteration:",l,"Model without contextual effects, Gy not observed, same draw and GX, GGX are instruments"))
    #print(lest1.2.1.2[l,])
    
    #different draw
    #We use GX constructed as instruments
    est1.2.2.1<-ivreg(Y1all ~ Xall + GY1call | Xall +  GXcall)
    sest1.2.2.1<-summary(est1.2.2.1, diagnostic=TRUE)
    lest1.2.2.1[l,] <- c(sest1.2.2.1$coefficients[,1],sest1.2.2.1$diagnostics[,3],cor((GY1all-GY1call),GXcall))
    #print(paste("Iteration:",l,"Model without contextual effects, Gy not observed same draw and GX is instrument"))
    #print(lest1.2.2.1[l,])
    
    #We use GX and GGX constructed as instruments
    est1.2.2.2<-ivreg(Y1all ~ Xall + GY1call | Xall +  GXcall + G2Xcall)
    sest1.2.2.2<-summary(est1.2.2.2, diagnostic=TRUE)
    lest1.2.2.2[l,] <- c(sest1.2.2.2$coefficients[,1],sest1.2.2.2$diagnostics[,3],cor((GY1all-GY1call),cbind(GXcall,G2Xcall)))
    #print(paste("Iteration:",l,"Model without contextual effects, Gy not observed, same draw and GX, GGX are instruments"))
    #print(lest1.2.2.2[l,])
    print(paste0("lambda: ", i,"/",length(lambda),  " - Iteration: ", l, "/", Iteration))
  }
  
  Lest1[[i]] <- list(lest1.1.1,lest1.1.2,lest1.2.1.1,lest1.2.1.2,lest1.2.2.1,lest1.2.2.2)
}

save(Lest1, lambda, file = "Results/Lest1.Rdata")
################################################
rm(list = ls())
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
load("Results/Lest1.Rdata")

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
temp <- Lest1[[which(round(lambda, 2) == lambdai)]]

lambda[which(round(lambda, 2) == lambdai)]

lest1.1.1 <- temp[[1]]
lest1.1.2 <- temp[[2]]
lest1.2.1.1 <- temp[[3]]
lest1.2.1.2 <- temp[[4]]
lest1.2.2.1 <- temp[[5]]
lest1.2.2.2 <- temp[[6]]

(sum.lest1.1.1 <- round(t(as.matrix(apply(lest1.1.1 , 2, my.sum))), 3))
(sum.lest1.1.2 <- round(t(as.matrix(apply(lest1.1.2 , 2, my.sum))), 3))
(sum.lest1.2.1.1 <- round(t(as.matrix(apply(lest1.2.1.1 , 2, my.sum))), 3))
(sum.lest1.2.1.2 <- round(t(as.matrix(apply(lest1.2.1.2 , 2, my.sum))), 3))
(sum.lest1.2.2.1 <- round(t(as.matrix(apply(lest1.2.2.1 , 2, my.sum))), 3))
(sum.lest1.2.2.2 <- round(t(as.matrix(apply(lest1.2.2.2 , 2, my.sum))), 3))

write.csv(sum.lest1.1.1, file = paste("Results/sum.lest1.1.1_", lambdai, ".csv"))
write.csv(sum.lest1.1.2, file = paste("Results/sum.lest1.1.2_", lambdai, ".csv"))
write.csv(sum.lest1.2.1.1, file = paste("Results/sum.lest1.2.1.1_", lambdai, ".csv"))
write.csv(sum.lest1.2.1.2, file = paste("Results/sum.lest1.2.1.2_", lambdai, ".csv"))
write.csv(sum.lest1.2.2.1, file = paste("Results/sum.lest1.2.2.1_", lambdai, ".csv"))
write.csv(sum.lest1.2.2.2, file = paste("Results/sum.lest1.2.2.2_", lambdai, ".csv"))


#### Graph
alpha.mean.1.1.1 <-
  unlist(lapply(lindex, function(w) {
    mean(Lest1[[w]][[1]][, "alpha"])
  }))
alpha.1st.Qu.1.1.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest1[[w]][[1]][, "alpha"], 0.25)
  }))
alpha.2st.Qu.1.1.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest1[[w]][[1]][, "alpha"], 0.75)
  }))
alpha.INF.CI.1.1.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest1[[w]][[1]][, "alpha"], 0.025)
  }))
alpha.SUP.CI.1.1.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest1[[w]][[1]][, "alpha"], 0.975)
  }))

tab.graph.1.1.1 <-
  data.frame(
    lambda = lambda[lindex],
    mean = alpha.mean.1.1.1,
    Inf.IQI = alpha.1st.Qu.1.1.1,
    Sup.IQI = alpha.2st.Qu.1.1.1,
    Inf.CI = alpha.INF.CI.1.1.1,
    Sup.CI = alpha.SUP.CI.1.1.1
  )

alpha.mean.1.2.2.1 <-
  unlist(lapply(lindex, function(w) {
    mean(Lest1[[w]][[5]][, "alpha"])
  }))
alpha.1st.Qu.1.2.2.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest1[[w]][[5]][, "alpha"], 0.25)
  }))
alpha.2st.Qu.1.2.2.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest1[[w]][[5]][, "alpha"], 0.75)
  }))
alpha.INF.CI.1.2.2.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest1[[w]][[5]][, "alpha"], 0.025)
  }))
alpha.SUP.CI.1.2.2.1 <-
  unlist(lapply(lindex, function(w) {
    quantile(Lest1[[w]][[5]][, "alpha"], 0.975)
  }))

tab.graph.1.2.2.1 <-
  data.frame(
    lambda = lambda[lindex],
    mean = alpha.mean.1.2.2.1,
    Inf.IQI = alpha.1st.Qu.1.2.2.1,
    Sup.IQI = alpha.2st.Qu.1.2.2.1,
    Inf.CI = alpha.INF.CI.1.2.2.1,
    Sup.CI = alpha.SUP.CI.1.2.2.1
  )



g1 <- ggplot(tab.graph.1.1.1, aes(x = lambda)) + theme_bw() +
  geom_ribbon(aes(ymin = Inf.CI, ymax = Sup.CI, fill = "95% CI"), alpha = 1) +
  geom_ribbon(aes(ymin = Inf.IQI, ymax = Sup.IQI, fill = "Interquartile Interval"), alpha = 1) +
  geom_line(aes(y = mean, colour = "Mean")) + ylab("") +
  scale_fill_manual(NULL, values = c("#ccccff", "#55ffff")) +
  scale_colour_manual(NULL, values = "#ff0000") +
  theme(
    legend.position = c(0.3, 0.885),
    legend.box = 'horizontal',
    axis.ticks = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 1)) +
  ylab("alpha") +
  ggtitle(paste(expression(Gy), "is oberved")) +
  scale_y_continuous(breaks = seq(0.35, 0.45, 0.02),
                     limits = c(0.35, 0.45)) +
  scale_x_continuous(breaks = seq(0,2,0.5), labels = as.expression(seq(0,2,0.5)), 
                     sec.axis = sec_axis(~ ., name = "Correlation",
                                         labels = c(1,0.394,0.173,0.092,0.056)))

g2 <-
  ggplot(tab.graph.1.2.2.1, aes(x = lambda)) +  theme_bw() +
  geom_ribbon(aes(ymin = Inf.CI, ymax = Sup.CI, fill = "95% CI"), alpha = 1) +
  geom_ribbon(aes(ymin = Inf.IQI, ymax = Sup.IQI, fill = "Interquartile Interval"), alpha = 1) +
  geom_line(aes(y = mean, colour = "Mean")) + ylab("") +
  scale_fill_manual(NULL, values = c("#ccccff", "#55ffff")) +
  scale_colour_manual(NULL, values = "#ff0000") +
  theme(
    legend.position = c(0.3, 0.885),
    legend.box = 'horizontal',
    axis.ticks = element_blank()
  ) +
  ylab("") +
  ggtitle(paste(expression(Gy), "is not oberved")) +
  scale_y_continuous(breaks = seq(0.35, 0.45, 0.02),
                     limits = c(0.35, 0.45)) +
  scale_x_continuous(breaks = seq(0,2,0.5), labels = as.expression(seq(0,2,0.5)), 
                     sec.axis = sec_axis(~ ., name = "Correlation",
                                         labels = c(1,0.394,0.173,0.092,0.056)))

multiplot(g1, g2, cols = 2)
# size 11x4.5 inch