rm(list = ls())
library(readstata13)
library(ggplot2)
library(PartialNetwork)
library(Matrix)
library(DataCombine)
library(stargazer)
library(dplyr)
setwd("PartialNetwork")

# the finale data set is save in the 'filname' path and can be loaded if saved before
# otherwise, the code will be ran to prepare the data
filname        <- "DataPartialNetwork.rda"  
if (file.exists(filname)) {
  load(file = filname)
} else {
  # Data
  mydata       <- read.dta13("cleandta.dta")  # data from Stata
  mydata       <- mydata[order(mydata$sschlcde),]
  
  mislist      <- c(55555555, 77777777, 88888888, 99999999, 99959995)
  mf_coln      <- paste0("mf", 1:5, "aid")
  ff_coln      <- paste0("ff", 1:5, "aid")
  f_coln       <- c(mf_coln, ff_coln)
  
  # mislist is an ID
  if (sum(mydata$aid %in% mislist) > 0) {
    stop("mislist is an ID")
  } else {
    cat("mislist is not an ID: OK", "\n")
  }
  
  # list of variable (excluding reference variables for identification)
  va.names      <- c("female", "hispanic", "raceblack",   "raceasian",   "raceother", "melhigh", "memhigh",
                     "memiss", "mjprof", "mjother",  "mjmiss", "age", "gpa")
  
  # Are there NAs?
  apply(mydata[,va.names], 2, function(w) sum(is.na(w)))
  
  # Keep row without NA
  keep1         <- as.logical((rowSums(is.na(mydata[,va.names])) == 0))
  mydata        <- mydata[keep1,]
  
  # remove friend from different groups
  # remove self friendship
  # remove friend non found
  N       <- nrow(mydata)
  dscf    <- rep(0,N)
  sfre    <- rep(0,N)
  nffr    <- rep(0,N)
  for (i in 1:N) {
    for (j in f_coln) {
      k  <- which(mydata$aid == mydata[i, j, drop = TRUE])
      # remove if different school
      if (length(k) != 0) {
        if(mydata[i, "sschlcde", drop = TRUE] != mydata[k, "sschlcde", drop = TRUE]) {
          mydata[i, j]   <- -1
          dscf[i]        <- dscf[i] + 1
        }
        # remove if self frienship
        if(mydata[i, "aid", drop = TRUE] == mydata[k, "aid", drop = TRUE]) {
          mydata[i, j]   <- -2
          sfre[i]        <- sfre[i] + 1
        }
      }
      else {
        if (!((mydata[i, j, drop = TRUE] %in% mislist) | is.na(mydata[i, j, drop = TRUE]))) {
          mydata[i, j]   <- -3
          nffr[i]        <- nffr[i] + 1
        }
      }
    }
  }
  
  cat("remove", sum(dscf), "link(s) because students from different schools: their code are recode as -1", "\n")
  cat("remove", sum(sfre), "self-friendship(s): their code are recoded as -2", "\n")
  cat("remove", sum(nffr), "non-found friends: their code are recoded as -3", "\n")
  rm(list = c("i", "j", "k"))
  
  # School size max
  sch.size.max  <- 200
  
  # Keep schools < School size max
  schtable      <- table(mydata$sschlcde)
  school        <- as.numeric(names(schtable))
  sch.size      <- as.numeric(schtable)
  school        <- school[sch.size < sch.size.max]
  sch.size      <- sch.size[sch.size < sch.size.max]
  nsch          <- length(school)
  
  # Update data
  mydata        <- mydata[mydata$sschlcde %in% school,]
  
  # dependent variable
  mydata$y     <- mydata[,tail(va.names,1), drop = TRUE]
  
  mislistmis   <- c(55555555, 99999999, 99959995)
  
  
  # This function prepares the data and the network 
  gen.data  <- function(db) {
    G       <- vector("list", nsch)
    Gobmis  <- vector("list", nsch)
    Gobtop  <- vector("list", nsch)
    Gobtmis <- vector("list", nsch)
    nmatchm <- vector("list", nsch)
    nmatchf <- vector("list", nsch)
    matchm  <- vector("list", nsch)
    matchf  <- vector("list", nsch)
    friendm <- vector("list", nsch)
    friendf <- vector("list", nsch)
    gendf   <- vector("list", nsch)
    X       <- as.matrix(db[, va.names[-length(va.names)]])    
    Y       <- db[,"y", drop = TRUE]
    
    
    # output for this loop
    # G, 
    # Gob (the observed part of G)
    # X containing also number of missing links
    # ni is number of students in all schools
    for (i in 1:nsch) {
      cat("School :", i, "/", nsch, "\n")
      schi         <- school[i]
      dbi          <- db %>% filter(sschlcde == schi)
      Ni           <- nrow(dbi)
      Gi           <- matrix(0, Ni, Ni)
      Giom         <- matrix(1, Ni, Ni)
      Giot         <- matrix(1, Ni, Ni)
      Giotm        <- matrix(1, Ni, Ni)
      gendfi       <- matrix(0, Ni, Ni)
      
      nmatchim     <- numeric() #will contain missing male links
      nmatchif     <- numeric() #will contain missing female links
      matchim      <- numeric() #will contain male links
      matchif      <- numeric() #will contain female links
      for (j in 1:Ni) {
        idxm       <- which(dbi$aid %in% unlist(dbi[j, mf_coln]))
        idxf       <- which(dbi$aid %in% unlist(dbi[j, ff_coln]))
        matchim[j] <- length(idxm) # observed male links
        matchif[j] <- length(idxf) # observed female links
        idx        <- c(idxm, idxf)
        Gi[j, idx] <- 1
        
        # missing links
        idxm.miss  <- which(unlist(dbi[j, mf_coln]) %in% c(mislistmis, -3))
        idxf.miss  <- which(unlist(dbi[j, ff_coln]) %in% c(mislistmis, -3))
        nmatchim[j]<- length(idxm.miss) # Number of missing male links
        nmatchif[j]<- length(idxf.miss) # Number of missing female links
        
        # If there are messing links, then we doubt the values
        if(nmatchim[j] > 0) {
          # Giom[j, (Gi[j,] == 0)] <- 0
          Giom[j, dbi$female == 0] <- 0
        }
        if(nmatchif[j] > 0) {
          # Giom[j, (Gi[j,] == 0)] <- 0
          Giom[j, dbi$female == 1] <- 0
        }        
        
        # Top coding
        # If top coding for any sex, then we doubt all zeros in Gi[j,] for the same sex
        # male
        if ((length(idxm) + nmatchim[j]) == 5){# There are 5 male friends declared
          Giot[j, (Gi[j,] == 0) & (dbi$female == 0)] <- 0
        }
        # female
        if ((length(idxf) + nmatchif[j]) == 5){# There are 5 female friends declared
          Giot[j, (Gi[j,] == 0) & (dbi$female == 1)] <- 0
        }
        
        # Missing and top coding
        Giotm[j,]  <- Giom[j,]*Giot[j,]
        
        #gender 
        gendfi[j,] <- dbi$female
      }
      
      #  store G
      diag(Gi)     <- 0
      diag(Giom)   <- 0
      diag(Giot)   <- 0
      diag(Giotm)  <- 0
      G[[i]]       <- Gi
      Gobmis[[i]]  <- Giom
      Gobtop[[i]]  <- Giot
      Gobtmis[[i]] <- Giotm
      
      # matched
      matchm[[i]]  <- matchim
      matchf[[i]]  <- matchif
      
      # unmatched
      nmatchm[[i]] <- nmatchim
      nmatchf[[i]] <- nmatchif
      
      # friends 
      friendm[[i]] <- matchim + nmatchim
      friendf[[i]] <- matchif + nmatchif
      
      # gender 
      gendf[[i]]   <- gendfi
    }
    
    list(G       = G,
         Gobmis  = Gobmis,
         Gobtop  = Gobtop,
         Gobtmis = Gobtmis,
         X       = X,
         matchm  = matchm,
         matchf  = matchf,
         nmatchm = nmatchm,
         nmatchf = nmatchf,
         friendm = friendm,
         friendf = friendf,
         Y       = Y,
         gendf   = gendf)
    
  }
  
  # Use the function to prepare the data
  tmp     <- gen.data(mydata)
  G       <- tmp$G
  Gobmis  <- tmp$Gobmis
  Gobtop  <- tmp$Gobtop
  Gobtmis <- tmp$Gobtmis
  X       <- tmp$X
  matchm  <- tmp$matchm
  matchf  <- tmp$matchf
  nmatchm <- tmp$nmatchm
  nmatchf <- tmp$nmatchf
  friendm <- tmp$friendm
  friendf <- tmp$friendf
  Y       <- tmp$Y
  gendf   <- tmp$gendf
  
  
  # dataset for the logistic model
  # variable to include 1 (indicator == 1 if same value)
  va.log1       <- c("female", "hispanic", "racewhite", "raceblack", "raceasian", "melhigh", "memhigh",
                     "mjprof")
  # variable to include  2 (absolute value of the difference)
  va.log2       <- c("age")
  # distance
  dist1         <- function(x, y) as.numeric(x == y)
  dist2         <- function(x, y) abs(x - y)
  # data
  tmp           <- c(0, cumsum(sch.size))
  # for va.log1
  X1tmp         <- do.call("cbind", lapply(va.log1, function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx       <- mydata[c(tmp[x] + 1):tmp[x+1], z, drop = TRUE]   
      matrix(kronecker(va.zx, va.zx, FUN = dist1), sch.size[x])}))
  }))
  # for va.log2
  X2tmp         <- do.call("cbind", lapply(va.log2, function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx     <- mydata[c(tmp[x] + 1):tmp[x+1], z, drop = TRUE]   
      matrix(kronecker(va.zx, va.zx, FUN = dist2), sch.size[x])}))
  }))
  Xlogit            <- as.data.frame(cbind(X1tmp, X2tmp)) 
  colnames(Xlogit)  <- c("same.sex", va.log1[-1], "diff.age")
  
  save(list = ls(all = TRUE), file = filname)
}
############################ Descriptive stat ##################################
# Descriptive statistic function 
my.stat.des <- function(x) {
  out <- c(mean(x), sd(x), quantile(x, probs = c(0.005, 0.995, 0.025, 0.975, 0.05, 0.95, 0.25, 0.75)))
  names(out) <- c("Mean", "St. Dev.",   "Pctl(0.5)",   "Pctl(99.5)",   "Pctl(2.5)",   "Pctl(97.5)",   "Pctl(5)",   "Pctl(95)",   "Pctl(25)",   "Pctl(75)")
  out
}

# list of variable (excluding reference variables)
va.all.names   <- c("female", "hispanic", "racewhite", "raceblack",   "raceasian",   "raceother", "mehigh",
                    "melhigh", "memhigh", "memiss", "mjhome",  "mjprof", "mjother",  "mjmiss", "age", "gpa")
# all the variables 

allVar        <- mydata[,va.all.names]

# Descriptive stat
# This reproduced Table 7
sdes           <- round(t(apply(allVar, 2, my.stat.des)), 3)[,c(1,2,9,10)]
sdes
print(sdes)

# convert result into latex
c1.names   <- c("Female", "Hispanic", "Race", "", "", "" , "", "Mother Education", "", "", "", "",
                "Mother Job", "", "", "", "", "Age", "GPA")
c2.names   <- c("", "", "", "White", "Black", "Asian", "Other", "", "high", "< high", "> high",
                "Missing", "", "Stay home",  "Professional", "Other", "Missing", "", "")

sdes       <- as.data.frame(cbind(format(round(sdes[,c(1,2)],3), 3), sdes[,-c(1,2)]))

for (i in 1:4) {
  sdes[,i] <- as.character(sdes[,i])
}

sdes <- InsertRow(sdes, NewRow = rep("NA", 4), RowNum = 3)
sdes <- InsertRow(sdes, NewRow = rep("NA", 4), RowNum = 8)
sdes <- InsertRow(sdes, NewRow = rep("NA", 4), RowNum = 13)

sdes <- as.data.frame(cbind(c1.names, c2.names, sdes))   

for (i in 3:6) {
  sdes[grepl("NA", sdes[,i]),i]     <- ""
}

stargazer(sdes, summary = FALSE, rownames = FALSE)

# number of observed friends
matchmall     <- unlist(matchm)
matchfall     <- unlist(matchf)
matchall      <- matchmall + matchfall 
summary(matchall - unlist(lapply(G, rowSums))) #should be zero
sum(matchall)

# number of unmatched friends
nmatchmall    <- unlist(nmatchm)
nmatchfall    <- unlist(nmatchf)
nmatchall     <- nmatchmall + nmatchfall 
sum(nmatchall)

# declared friends
friendmall    <- unlist(friendm)
friendfall    <- unlist(friendf)
friendall     <- friendmall + friendfall
sum(friendall)
summary(friendmall - matchmall - nmatchmall) #should be zero
summary(friendfall - matchfall - nmatchfall) #should be zero

# graph missing links
# This reproduces the Figure 4
ggplot(data = data.frame(mm = nmatchall), aes(x = mm)) +
  geom_bar(color = "black", fill = "#eeeeee") + 
  theme_bw() + xlab("") + ylab("Frequency") + 
  scale_x_discrete(limits = 0:10) 
# size 7 × 4 inch

# Proportion of observed 
sum(sapply(Gobmis, function(x) sum(x == 1)))/sum(sch.size*(sch.size - 1))
# proportion of missing to be inferred
sum(sapply(Gobmis, function(x) sum(x == 0) - nrow(x)))/sum(sch.size*(sch.size - 1))
# Proportion of missing to be inferred due to error code and top coding
sum(sapply(Gobtmis, function(x) sum(x == 0) - nrow(x)))/sum(sch.size*(sch.size - 1))

################################ Estimation #############################
######## dataset and model settings
dataset    <- data.frame(GPA = Y, X)
Xlogit     <- as.data.frame(Xlogit)
exp.var    <- c("female", "hispanic", "raceblack", "raceasian", "raceother", "melhigh", "memhigh", 
                "memiss", "mjprof", "mjother", "mjmiss", "age")
Model      <- as.formula(paste0("GPA ~ ", paste0(exp.var, collapse = "+")))
Kv         <- 2*ncol(X) + 1  # Number of exogenous explanatory variables

# Hyperparameters
hyperp     <- list("mutheta"       = rep(0,Kv),
                   "invstheta"     = diag(Kv)/100,
                   "muzeta"        = 0,
                   "invszeta"      = 2,
                   "a"             = 4.2,
                   "b"             = (4.2 - 2)*0.5)
####################################  Observed network as given
# Beyesian estimator
set.seed(123)
obs.bayes.est    <- mcmcSAR(formula    = Model,
                            contextual = TRUE,
                            G0.obs     = "all",
                            G0         = G,
                            hyperparms = hyperp,
                            data       = dataset,
                            iteration  = 2e4)
saveRDS(obs.bayes.est, file = "obs.bayes.est.RDS")
summary(obs.bayes.est)
plot(obs.bayes.est, mar = c(2, 2, 1, 1))
plot(obs.bayes.est, plot.type = "dens", mar = c(2, 2, 1, 1))

# write.csv(summary(obs.bayes.est)$posterior$theta, file = "tmp.csv")

# SGMM estimator
# W = I
set.seed(123)
obs.sgmm.est1    <- smmSAR(formula    = Model,
                           contextual = TRUE,
                           dnetwork   = G, # there is no missing link 
                           smm.ctr    = list(R = 1L, iv.power = 2L, opt.tol = 1e-7, print = TRUE),
                           data       = dataset)
summary(obs.sgmm.est1)
saveRDS(obs.sgmm.est1, file = "obs.sgmm.est1.RDS")

# write.csv(print(summary(obs.sgmm.est1))$coefficients, file = "tmp.csv")

# W = (Z'Z)^{-1}
Gsimnorm         <- norm.network(G)
GXsim            <- peer.avg(Gsimnorm, dataset[,exp.var])
GGXsim           <- peer.avg(Gsimnorm, GXsim)
W                <- solve(crossprod(as.matrix(cbind(1, dataset[,exp.var], GXsim, GGXsim))))
set.seed(123)
obs.sgmm.est2    <- smmSAR(formula    = Model,
                           contextual = TRUE,
                           dnetwork   = G, # there is no missing link 
                           W          = W,
                           smm.ctr    = list(R = 1L, iv.power = 2L, opt.tol = 1e-7, print = TRUE),
                           data       = dataset)
summary(obs.sgmm.est2)
saveRDS(obs.sgmm.est2, file = "obs.sgmm.est2.RDS")
# write.csv(print(summary(obs.sgmm.est2))$coefficients, file = "tmp.csv")

####################################  Reconstructed network
######################### only missing links, we ignore top coding here
# Weights to fix the selection issue
tabfriends <- data.frame(totfriendm = matchmall + nmatchmall, 
                         totfriendf = matchfall + nmatchfall,
                         misfriendm = nmatchmall,
                         misfriendf = nmatchfall)

# For each student, we only consider their males (resp. females) friends if there are no missing males (resp. females) links
# This is who we compute the weigh to address the selection issue.
weightm <- unlist(tabfriends %>% group_by(totfriendm) %>% 
                    summarise(w = length(misfriendm)/sum(misfriendm == 0)) %>% select(w))
weightf <- unlist(tabfriends %>% group_by(totfriendf) %>% 
                    summarise(w = length(misfriendf)/sum(misfriendf == 0)) %>% select(w))
weights <- lapply(1:nsch, function(x) {
  wm    <- matrix(rep(weightm[friendm[[x]] + 1], each = sch.size[x]), sch.size[x], byrow = TRUE)
  wf    <- matrix(rep(weightf[friendf[[x]] + 1], each = sch.size[x]), sch.size[x], byrow = TRUE)
  wf*gendf[[x]] + wm*(1 - gendf[[x]])
})
weights <- mat.to.vec(weights)[as.logical(mat.to.vec(Gobmis))]

# Beyesian estimator
set.seed(123)
mlinks           <- list(model          = "logit", 
                         mlinks.formula = paste0("~ ", paste0(c("same.sex", va.log1[-1], "diff.age"), collapse = "+")), 
                         mlinks.data    = Xlogit,
                         weights        = weights)

miss.bayes.est   <- mcmcSAR(formula    = Model,
                            contextual = TRUE,
                            G0.obs     = Gobmis,
                            G          = G,
                            hyperparms = hyperp,
                            ctrl.mcmc  = list(print.level = 2),
                            mlinks     = mlinks,
                            data       = dataset,
                            iteration  = 2e4)
saveRDS(miss.bayes.est, file = "miss.bayes.est.RDS")
summary(miss.bayes.est)
plot(miss.bayes.est, mar = c(2, 2, 1, 1))

# write.csv(rbind(summary(miss.bayes.est)$posterior$theta, summary(miss.bayes.est)$posterior$rho), file = "tmp.csv")
# SGMM estimator
# We first estimate the distribution of the network using a logit regression on the observed network
Aobs.mis          <- as.logical(mat.to.vec(lapply(G, function(x) 1*(x > 0))))[as.logical(mat.to.vec(Gobmis))]
Xlogit.mis        <- as.matrix(Xlogit[as.logical(mat.to.vec(Gobmis)),])
miss.logit        <- glm(Aobs.mis ~ Xlogit.mis, family = binomial(link = "logit"), weights = weights)

summary(miss.logit)
rho.mis           <- miss.logit$coefficients
var.rho           <- summary(miss.logit)$cov.unscaled
saveRDS(miss.logit, file = "miss.logit.RDS")
# write.csv(summary(miss.logit)$coefficients, file = "tmp.csv")

# estimated distribution for unobserved links
dnetwork.miss     <- vec.to.mat(plogis(c(cbind(1, as.matrix(Xlogit))%*%rho.mis)), N = sch.size)
dnetwork.miss     <- lapply(1:nsch, function(x) (G[[x]] > 0)*Gobmis[[x]] + dnetwork.miss[[x]]*(1 - Gobmis[[x]]))

# W = I
set.seed(123)
miss.sgmm.est1     <- smmSAR(formula    = Model,
                             contextual = TRUE,
                             dnetwork   = dnetwork.miss, 
                             smm.ctr    = list(R = 1000L, iv.power = 2L, opt.tol = 1e-4, print = TRUE),
                             data       = dataset)
summary(miss.sgmm.est1)
saveRDS(miss.sgmm.est1, file = "miss.sgmm.est1.RDS")

# variance computation
fdist.miss        <- function(rho.mis, var.rho, Xlogit, G, Gobmis){
  M               <- length(G)
  N               <- sapply(G, nrow)
  rho.sim         <- MASS::mvrnorm(mu = rho.mis, Sigma = var.rho)
  dsim            <- vec.to.mat(plogis(c(cbind(1, as.matrix(Xlogit))%*%rho.sim)), N = N)
  lapply(1:M, function(x) (G[[x]] > 0)*Gobmis[[x]] + dsim[[x]]*(1 - Gobmis[[x]]))
}
fdist_args        <- list(rho.mis = rho.mis, var.rho = var.rho, Xlogit = Xlogit, G = G, Gobmis = Gobmis)
smiss.sgmm.est1    <- summary(miss.sgmm.est1, dnetwork = dnetwork.miss, data = dataset, 
                              .fun = fdist.miss, .args = fdist_args, sim = 500L, ncores = 5L)
print(smiss.sgmm.est1)
saveRDS(smiss.sgmm.est1, file = "smiss.sgmm.est1.RDS")

# write.csv(print(smiss.sgmm.est1)$coefficients, file = "tmp.csv")

# W = (Z'Z)^{-1}
Gsimnorm         <- norm.network(sim.network(dnetwork.miss))
GXsim            <- peer.avg(Gsimnorm, dataset[,exp.var])
GGXsim           <- peer.avg(Gsimnorm, GXsim)
W                <- solve(crossprod(as.matrix(cbind(1, dataset[,exp.var], GXsim, GGXsim))))
set.seed(123)
miss.sgmm.est2     <- smmSAR(formula    = Model,
                             contextual = TRUE,
                             dnetwork   = dnetwork.miss, 
                             W          = W,
                             smm.ctr    = list(R = 1000L, iv.power = 2L, opt.tol = 1e-4, print = TRUE),
                             data       = dataset)
saveRDS(miss.sgmm.est2, file = "miss.sgmm.est2.RDS")
summary(miss.sgmm.est2)
smiss.sgmm.est2    <- summary(miss.sgmm.est2, dnetwork = dnetwork.miss, data = dataset, 
                              .fun = fdist.miss, .args = fdist_args, sim = 500L, ncores = 5L)
print(smiss.sgmm.est2)
saveRDS(smiss.sgmm.est2, file = "smiss.sgmm.est2.RDS")

# write.csv(print(smiss.sgmm.est2)$coefficients, file = "tmp.csv")

######################### Missing links and top coding
# We first fit the number of friends
censure           <- (friendmall >= 5) | (friendfall >= 5)
Xpoisson          <- as.matrix(cbind(fastDummies::dummy_cols(mydata$sschlcde)[,-1], 
                                     mydata[,c("female", "hispanic", "racewhite", "raceblack", "raceasian", 
                                               "melhigh", "memhigh", "mjprof", "age")]))
lcensure          <- friendall
rcensure          <- ifelse(friendmall >= 5, unlist(lapply(gendf, function(x) rowSums(x == 0))), friendmall)
rcensure          <- rcensure + ifelse(friendfall >= 5, unlist(lapply(gendf, function(x) rowSums(x == 1))), friendfall)

par               <- runif(ncol(Xpoisson), -1, 1)
poisson.optin1    <- optim(par, f_rcpoisson, control = list(maxit = 1e9, reltol = 1e-20), y = friendall, 
                           X = Xpoisson, censure = censure, lcensure = lcensure, rcensure = rcensure)

poisson.optin2    <- nlm(p = poisson.optin1$par, f = f_rcpoisson, iterlim = 1e9, gradtol = 1e-20, y = friendall, 
                         X = Xpoisson, censure = censure, lcensure = lcensure, rcensure = rcensure)

est.nfriends      <- expy(beta = poisson.optin2$estimate, X = Xpoisson, lcensure = lcensure, rcensure = rcensure)
est.nfriendm      <- ifelse(friendmall < 5, friendmall, ifelse(friendfall < 5, est.nfriends - friendfall, est.nfriends/2))
est.nfriendf      <- ifelse(friendfall < 5, friendfall, ifelse(friendmall < 5, est.nfriends - friendmall, est.nfriends/2))
save(est.nfriends, est.nfriendm, est.nfriendf, file = "~/Dropbox/Papers - In progress/Partial Network/AddHealth/est.nfriend.rda")

# Weights to fix the selection issue
# We use all the student. Links are weighed to address the selection issue
load("est.nfriend.rda")
est.nfriendm <- floor(est.nfriendm)
est.nfriendf <- floor(est.nfriendf)

# We first compute weights for missing links and we correct for top coding
# Weights to fix the selection issue
tabfriends <- data.frame(totfriendm = matchmall + nmatchmall, 
                         totfriendf = matchfall + nmatchfall,
                         misfriendm = nmatchmall,
                         misfriendf = nmatchfall)

# For each student, we only consider their males (resp. females) friends if there are no missing males (resp. females) links
# This is who we compute the weigh to address the selection issue.
weightm <- unlist(tabfriends %>% group_by(totfriendm) %>% 
                    summarise(w = length(misfriendm)/sum(misfriendm == 0)) %>% select(w))
weightf <- unlist(tabfriends %>% group_by(totfriendf) %>% 
                    summarise(w = length(misfriendf)/sum(misfriendf == 0)) %>% select(w))
weights <- lapply(1:nsch, function(x) {
  wm    <- matrix(rep(weightm[friendm[[x]] + 1], each = sch.size[x]), sch.size[x], byrow = TRUE)
  wf    <- matrix(rep(weightf[friendf[[x]] + 1], each = sch.size[x]), sch.size[x], byrow = TRUE)
  # we correct for top coding
  nmx   <- sum(gendf[[x]][1,] == 0)
  nfx   <- sum(gendf[[x]][1,] == 1)
  csum  <- c(0, cumsum(sch.size))
  idx   <- (csum[x] + 1):csum[x + 1]
  tcm1  <- matrix(rep(est.nfriendm[idx]/friendm[[x]], each = sch.size[x]), sch.size[x], byrow = TRUE)
  tcm0  <- matrix(rep((nmx - 1 - est.nfriendm[idx])/(nmx - 1 - friendm[[x]]), each = sch.size[x]), sch.size[x], byrow = TRUE)
  tcf1  <- matrix(rep(est.nfriendf[idx]/friendf[[x]], each = sch.size[x]), sch.size[x], byrow = TRUE)
  tcf0  <- matrix(rep((nfx - 1 - est.nfriendf[idx])/(nfx - 1 - friendf[[x]]), each = sch.size[x]), sch.size[x], byrow = TRUE)
  tcm1[is.na(tcm1)] <- 0
  tcm0[is.na(tcm0)] <- 0
  tcf1[is.na(tcf1)] <- 0
  tcf0[is.na(tcf0)] <- 0
  tcm   <- tcm1*G[[x]] + tcm0*(1 - G[[x]])
  tcf   <- tcf1*G[[x]] + tcf0*(1 - G[[x]])
  tcf*wf*gendf[[x]] + tcm*wm*(1 - gendf[[x]])
})
weights <- mat.to.vec(weights)[as.logical(mat.to.vec(Gobmis))]

# Network distribution
Aobs.tmis         <- as.logical(mat.to.vec(lapply(G, function(x) 1*(x > 0))))[as.logical(mat.to.vec(Gobmis))]
Xlogit.tmis       <- as.matrix(Xlogit[as.logical(mat.to.vec(Gobmis)),])
tmiss.logit       <- glm(Aobs.tmis ~ Xlogit.tmis, family = binomial(link = "logit"), weights = weights)

summary(tmiss.logit)
rho.tmis          <- tmiss.logit$coefficients
var.rho           <- summary(tmiss.logit)$cov.unscaled
saveRDS(tmiss.logit, file = "tmiss.logit.RDS")
# write.csv(summary(tmiss.logit)$coefficients, file = "tmp.csv")

# Beyesian estimator
set.seed(123)
mlinks           <- list(model          = "logit", 
                         mlinks.formula = paste0("~ ", paste0(c("same.sex", va.log1[-1], "diff.age"), collapse = "+")), 
                         mlinks.data    = Xlogit,
                         estimates      = list(mu = rho.tmis, var.rho = var.rho))

tmiss.bayes.est  <- mcmcSAR(formula    = Model,
                            contextual = TRUE,
                            G0.obs     = Gobtmis,
                            G          = G,
                            hyperparms = hyperp,
                            ctrl.mcmc  = list(print.level = 2),
                            mlinks     = mlinks,
                            data       = dataset,
                            iteration  = 2e4)
saveRDS(tmiss.bayes.est, file = "tmiss.bayes.est.RDS")
summary(tmiss.bayes.est)
plot(tmiss.bayes.est, mar = c(2, 2, 1, 1))

# write.csv(rbind(summary(tmiss.bayes.est)$posterior$theta, summary(tmiss.bayes.est)$posterior$rho), file = "tmp.csv")

# SGMM estimator
dnetwork.tmiss    <- vec.to.mat(plogis(c(cbind(1, as.matrix(Xlogit))%*%rho.tmis)), N = sch.size)
dnetwork.tmiss    <- lapply(1:nsch, function(x) (G[[x]] > 0)*Gobtmis[[x]] + dnetwork.tmiss[[x]]*(1 - Gobtmis[[x]]))

# W = I
set.seed(123)
tmiss.sgmm.est1   <- smmSAR(formula    = Model,
                            contextual = TRUE,
                            dnetwork   = dnetwork.tmiss, 
                            smm.ctr    = list(R = 1000L, iv.power = 2L, opt.tol = 1e-4, print = TRUE),
                            data       = dataset)
summary(tmiss.sgmm.est1)
saveRDS(tmiss.sgmm.est1, file = "tmiss.sgmm.est1.RDS")

# variance computation
fdist.tmiss       <- function(rho.tmis, var.rho, Xlogit, G, Gobtmis){
  M               <- length(G)
  N               <- sapply(G, nrow)
  rho.sim         <- MASS::mvrnorm(mu = rho.tmis, Sigma = var.rho)
  dsim            <- vec.to.mat(plogis(c(cbind(1, as.matrix(Xlogit))%*%rho.sim)), N = N)
  lapply(1:M, function(x) (G[[x]] > 0)*Gobtmis[[x]] + dsim[[x]]*(1 - Gobtmis[[x]]))
}
fdist_args        <- list(rho.tmis = rho.tmis, var.rho = var.rho, Xlogit = Xlogit, G = G, Gobtmis = Gobtmis)

stmiss.sgmm.est1  <- summary(tmiss.sgmm.est1, dnetwork = dnetwork.tmiss, data = dataset, 
                             .fun = fdist.tmiss, .args = fdist_args, sim = 500L, ncores = 5L)
saveRDS(stmiss.sgmm.est1, file = "stmiss.sgmm.est1.RDS")
print(stmiss.sgmm.est1)

# write.csv(print(stmiss.sgmm.est1)$coefficients, file = "tmp.csv")

# W = (Z'Z)^{-1}
Gsimnorm          <- norm.network(sim.network(dnetwork.tmiss))
GXsim             <- peer.avg(Gsimnorm, dataset[,exp.var])
GGXsim            <- peer.avg(Gsimnorm, GXsim)
W                 <- solve(crossprod(as.matrix(cbind(1, dataset[,exp.var], GXsim, GGXsim))))
set.seed(123)
tmiss.sgmm.est2   <- smmSAR(formula    = Model,
                            contextual = TRUE,
                            dnetwork   = dnetwork.tmiss, 
                            W          = W,
                            smm.ctr    = list(R = 1000L, iv.power = 2L, opt.tol = 1e-4, print = TRUE),
                            data       = dataset)
saveRDS(tmiss.sgmm.est2, file = "tmiss.sgmm.est2.RDS")
summary(tmiss.sgmm.est2)

stmiss.sgmm.est2  <- summary(tmiss.sgmm.est2, dnetwork = dnetwork.tmiss, data = dataset, 
                             .fun = fdist.tmiss, .args = fdist_args, sim = 500L, ncores = 5L)
saveRDS(stmiss.sgmm.est2, file = "stmiss.sgmm.est2.RDS")
print(stmiss.sgmm.est2)

# write.csv(print(stmiss.sgmm.est2)$coefficients, file = "tmp.csv")

####### Plot simulations
obs.ntw           <- obs.bayes.est$posterior
recons.ntw        <- tmiss.bayes.est$posterior$theta
form.ntw          <- tmiss.bayes.est$posterior$rho
Xnames            <- c("Female", "Hispanic", paste("Race =", c("Black", "Asian", "Other")),
                       paste0("Mother Edu ", c("< High", "> High", "= Missing")), 
                       paste("Mother Job =", c("Professional", "Other", "Missing")), "Age")

XnamesRho         <- c("Intercept", "Same sex", "Both Hispanic", "Both White", "Both Black",
                       "Both Asian", "Mums Educ < high", "Mums Educ > high",
                       "Mums Job Professional", "Age absolute diff")

### peer effect model, simulations
library(scales)
c1 = alpha("red", .8)
c2 = alpha("blue", .6)
par(fig = c(0, 1/6, 1 - 0.92/5, 1), mar = c(2, 2, 2, 2.1))
plot(obs.ntw[,26], type = "l", main = "Peer Effects", col = c1, xlab = "", ylab = "", ylim = c(min(obs.ntw[,26], recons.ntw[,26]), max(obs.ntw[,26], recons.ntw[,26])))
lines(recons.ntw[,26], type = "l", main = "Peer Effects", col = c2, xlab = "", ylab = "")

par(fig = c(1/6, 2/6, 1 - 0.92/5, 1), new = TRUE)
plot(obs.ntw[,1], type = "l", main = "Intercept", col = c1, xlab = "", ylab = "", ylim = c(min(obs.ntw[,1], recons.ntw[,1]), max(obs.ntw[,1], recons.ntw[,1])))
lines(recons.ntw[,1], type = "l", main = "Intercept", col = c2, xlab = "", ylab = "")

par(fig = c(2/6, 3/6, 1 - 0.92/5, 1), new = TRUE)
plot(obs.ntw[,27], type = "l", main = expression(sigma^2), col = c1, xlab = "", ylab = "", ylim = c(min(obs.ntw[,27], recons.ntw[,27]), max(obs.ntw[,27], recons.ntw[,27])))
lines(recons.ntw[,27], type = "l", main = expression(sigma^2), col = c2, xlab = "", ylab = "")

for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (1 + ceiling(i/6))*0.92/5 - 0.04, 1 - ceiling(i/6)*0.92/5 - 0.04),
      new = TRUE)
  plot(obs.ntw[,i + 1], type = "l", main = Xnames[i], col = c1, xlab = "", ylab = "", ylim = c(min(obs.ntw[,i + 1], recons.ntw[,i + 1]), max(obs.ntw[,i + 1], recons.ntw[,i + 1])))
  lines(recons.ntw[,i + 1], type = "l", main = Xnames[i], col = c2, xlab = "", ylab = "")
}

for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (3 + ceiling(i/6))*0.92/5 - 0.07999999, 1 - (2 + ceiling(i/6))*0.92/5 - 0.08),
      new = TRUE)
  plot(obs.ntw[,i + 13], type = "l", main = Xnames[i], col = c1, xlab = "", ylab = "", ylim = c(min(obs.ntw[,i + 13], recons.ntw[,i + 13]), max(obs.ntw[,i + 13], recons.ntw[,i + 13])))
  lines(recons.ntw[,i + 13], type = "l", main = Xnames[i], col = c2, xlab = "", ylab = "")
}

par(fig = c(0, 1, 1 - 2*0.92/5 - 0.04, 1 - 0.92/5), new = TRUE)
plot(x = 0, main = "Own Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

par(fig = c(0, 1, 1 - 4*0.92/5 - 0.08, 1 - 3*0.92/5 - 0.04), new = TRUE)
plot(x = 0, main = "Contextual Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

par(fig = c(0.5, 1, 0.75, 1), new = TRUE)
plot(x = 0, y = 0, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
legend(0, 1.2, legend=c("Observed Network", "Reconstructed Network"),
       col=c(c1, c2), lty = c(1,1), cex=1, box.lty=0)
 # 10x17

# peer effect model, densities
c1 = alpha("red", 1)
c2 = alpha("blue", 1)

par(fig = c(0, 1/6, 1 - 0.92/5, 1), mar = c(2, 2, 2, 2.1))
tmpo <- density(obs.ntw[2001:20000,26])
tmpr <- density(recons.ntw[2001:20000,26])
plot(tmpo, type = "l", main = "Peer Effects",
     col = c1, xlab = "", ylab = "", xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))))
lines(tmpr, col = c2, lty=2)

par(fig = c(1/6, 2/6, 1 - 0.92/5, 1), new = TRUE)
tmpo <- density(obs.ntw[2001:20000,1])
tmpr <- density(recons.ntw[2001:20000,1])
plot(tmpo, type = "l", main = "Intercept", col = c1, xlab = "", ylab = "", 
     xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))),
     ylim = c(min(c(tmpo$y, tmpr$y)), max(c(tmpo$y, tmpr$y))))
lines(tmpr, col = c2, lty=2)

par(fig = c(2/6, 3/6, 1 - 0.92/5, 1), new = TRUE)
tmpo <- density(obs.ntw[2001:20000,27])
tmpr <- density(recons.ntw[2001:20000,27])
plot(tmpo, type = "l", main = expression(sigma^2), col = c1, xlab = "", ylab = "",
     xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))),
     ylim = c(min(c(tmpo$y, tmpr$y)), max(c(tmpo$y, tmpr$y))))
lines(tmpr, col = c2, lty=2)

for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (1 + ceiling(i/6))*0.92/5 - 0.04, 1 - ceiling(i/6)*0.92/5 - 0.04),
      new = TRUE)
  tmpo <- density(obs.ntw[2001:20000, i + 1])
  tmpr <- density(recons.ntw[2001:20000, i + 1])
  plot(tmpo, type = "l", main = Xnames[i], col = c1, xlab = "", ylab = "",
       xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))),
       ylim = c(min(c(tmpo$y, tmpr$y)), max(c(tmpo$y, tmpr$y))))
  lines(tmpr, col = c2, lty=2)
}

for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (3 + ceiling(i/6))*0.92/5 - 0.07999999, 1 - (2 + ceiling(i/6))*0.92/5 - 0.08),
      new = TRUE)
  tmpo <- density(obs.ntw[2001:20000, i + 13])
  tmpr <- density(recons.ntw[2001:20000, i + 13])
  plot(tmpo, type = "l", main = Xnames[i], col = c1, xlab = "", ylab = "",
       xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))),
       ylim = c(min(c(tmpo$y, tmpr$y)), max(c(tmpo$y, tmpr$y))))
  lines(tmpr, col = c2, lty=2)
}
par(fig = c(0.5, 1, 0.75, 1), new = TRUE)
plot(x = 0, y = 0, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
legend(0, 1.2, legend=c("Observed Network", "Reconstructed Network"),
       col=c(c1, c2), lty=1:2, cex=1, box.lty=0)

par(fig = c(0, 1, 1 - 2*0.92/5 - 0.04, 1 - 0.92/5), new = TRUE)
plot(x = 0, main = "Own Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

par(fig = c(0, 1, 1 - 4*0.92/5 - 0.08, 1 - 3*0.92/5 - 0.04), new = TRUE)
plot(x = 0, main = "Contextual Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
# 10x17
### Network formation model
c1 = "blue"
par(mfrow = c(2, 5), mar = c(2, 2, 2, 2.1))
for (i in 1:10) {
  plot(form.ntw[,i], type = "l", main = XnamesRho[i], col = c1, xlab = "", ylab = "")
}
# 5x12

###################### Average number of friends
fnfriends   <- function(p.log, Gobs){
  network   <- 1*((as.matrix(cbind(1, Xlogit)) %*% p.log + rlogis(nrow(Xlogit))) > 0)
  network.n <- vec.to.mat(network, sch.size)
  network.n <- lapply(1:nsch, function(x) Gobs[[x]]*G[[x]] + (1 - Gobs[[x]])*network.n[[x]])
  unlist(lapply(network.n, rowSums))
}

library(doParallel)
summary(unlist(lapply(G, rowSums)))

# Missing
set.seed(123)
summary(apply(do.call(cbind, mclapply(10001:20000, function(x) fnfriends(miss.bayes.est$posterior$rho[x,], Gobmis), mc.cores = 20L)), 1, mean))
summary(fnfriends(miss.logit$coefficients, Gobmis))

# Missing and top coding
set.seed(123)
summary(apply(do.call(cbind, mclapply(10001:20000, function(x) fnfriends(tmiss.bayes.est$posterior$rho[x,], Gobtmis), mc.cores = 20L)), 1, mean))
summary(fnfriends(tmiss.logit$coefficients, Gobtmis))

# Peer effects
sobs.sgmm.est1 <- summary(obs.sgmm.est1)
dataint        <- data.frame(estimate = c(mean(obs.bayes.est$posterior[10001:20000, "Gy"]),
                                          mean(miss.bayes.est$posterior$theta[10001:20000, "Gy"]),
                                          mean(tmiss.bayes.est$posterior$theta[10001:20000, "Gy"]),
                                          sobs.sgmm.est1$estimates["Gy"],
                                          smiss.sgmm.est1$estimates["Gy"],
                                          stmiss.sgmm.est1$estimates["Gy"]),
                             LB       = c(quantile(obs.bayes.est$posterior[10001:20000, "Gy"], 0.025),
                                          quantile(miss.bayes.est$posterior$theta[10001:20000, "Gy"], 0.025),
                                          quantile(tmiss.bayes.est$posterior$theta[10001:20000, "Gy"], 0.025),
                                          sobs.sgmm.est1$estimates["Gy"] - 1.96*sqrt(sobs.sgmm.est1$cov[1,1]),
                                          smiss.sgmm.est1$estimates["Gy"] - 1.96*sqrt(smiss.sgmm.est1$cov[1,1]),
                                          stmiss.sgmm.est1$estimates["Gy"] - 1.96*sqrt(stmiss.sgmm.est1$cov[1,1])),
                             UB       = c(quantile(obs.bayes.est$posterior[10001:20000, "Gy"], 0.975),
                                          quantile(miss.bayes.est$posterior$theta[10001:20000, "Gy"], 0.975),
                                          quantile(tmiss.bayes.est$posterior$theta[10001:20000, "Gy"], 0.975),
                                          sobs.sgmm.est1$estimates["Gy"] + 1.96*sqrt(sobs.sgmm.est1$cov[1,1]),
                                          smiss.sgmm.est1$estimates["Gy"] + 1.96*sqrt(smiss.sgmm.est1$cov[1,1]),
                                          stmiss.sgmm.est1$estimates["Gy"] + 1.96*sqrt(stmiss.sgmm.est1$cov[1,1])),
                             Model    =  factor(c("Obsv.Bayes", "Miss.Bayes", "TopMiss.Bayes", "Obsv.SGMM", "Miss.SGMM", "TopMiss.SGMM"),
                                                levels = c("TopMiss.SGMM", "TopMiss.Bayes", "Miss.SGMM", "Miss.Bayes", "Obsv.SGMM", "Obsv.Bayes")))


ggplot(data = dataint, aes(y = Model, x = estimate)) + geom_errorbarh(aes(xmin = LB, xmax = UB, height = .3)) +
  geom_point(shape = 21, size = 2, fill = "white") +
  xlab("Peer effect estimate") + ylab("Model") + theme_bw()
#3*7
###################### Centrality
library(doParallel)
library(ggplot2)
# observed network
fcent.obs   <- function(alpha){
  unlist(lapply(norm.network(G), function(x) solve(diag(nrow(x)) - alpha*x, rep(1, nrow(x)))))
}
set.seed(123)
cent.obs    <- mclapply(10001:20000, function(x) fcent.obs(obs.bayes.est$posterior[x, "Gy"]), mc.cores = 20L)
scent.obs   <- apply(as.data.frame(cent.obs), 1, mean)

# reconstructed network
fcent.rec   <- function(p.log, alpha, Gobs){
  network   <- 1*((as.matrix(cbind(1, Xlogit)) %*% p.log + rlogis(nrow(Xlogit))) > 0)
  network.n <- vec.to.mat(network, sch.size)
  network.n <- lapply(1:nsch, function(x) Gobs[[x]]*G[[x]] + (1 - Gobs[[x]])*network.n[[x]])
  unlist(lapply(norm.network(network.n), function(x) solve(diag(nrow(x)) - alpha*x, rep(1, nrow(x)))))
}
set.seed(123)
cent.miss   <- mclapply(10001:20000, function(x) fcent.rec(miss.bayes.est$posterior$rho[x,], miss.bayes.est$posterior$theta[x, "Gy"], Gobmis), mc.cores = 20L)
set.seed(123)
cent.tmiss  <- mclapply(10001:20000, function(x) fcent.rec(tmiss.bayes.est$posterior$rho[x,], tmiss.bayes.est$posterior$theta[x, "Gy"], Gobtmis), mc.cores = 20L)
scent.miss  <- apply(as.data.frame(cent.miss), 1, mean)
scent.tmiss <- apply(as.data.frame(cent.tmiss), 1, mean)

centrality <- data.frame(observed = scent.obs, missing = scent.miss, top.missing = scent.tmiss)

saveRDS(centrality, file = "centrality.RDS")

# scatter plot
ggplot(data = centrality, aes(x = observed, y = top.missing)) + geom_point(colour = "#505050", size = 1) + theme_bw() + 
  xlab("Centrality: observed network") + ylab("Centrality: reconstructed network")

# boxplot
ggplot(data = centrality %>% filter(observed < 1.36), aes(y = factor(observed, labels = c("1.00", "1.35")), x = top.missing)) + 
  geom_boxplot() + ylab("Centrality: observed network") + xlab("Centrality: reconstructed network") + theme_bw()
#7*3
