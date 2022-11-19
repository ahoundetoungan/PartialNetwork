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
      k  <- which(mydata$aid == mydata[i, j])
      # remove if different school
      if (length(k) != 0) {
        if(mydata[i, "sschlcde"] != mydata[k, "sschlcde"]) {
          mydata[i, j]   <- -1
          dscf[i]        <- dscf[i] + 1
        }
        # remove if self frienship
        if(mydata[i, "aid"] == mydata[k, "aid"]) {
          mydata[i, j]   <- -2
          sfre[i]        <- sfre[i] + 1
        }
      }
      else {
        if (!((mydata[i, j] %in% mislist) | is.na(mydata[i, j]))) {
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
  sch.size.max <- 200
  
  # Keep schools < School size max
  
  schtable     <- table(mydata$sschlcde)
  school       <- as.numeric(names(schtable))
  sch.size     <- as.numeric(schtable)
  school       <- school[sch.size < sch.size.max]
  sch.size     <- sch.size[sch.size < sch.size.max]
  
  # Update data
  mydata       <- mydata[mydata$sschlcde %in% school,]
  nsch         <- length(school)
  
  # dependent variable
  mydata$y     <- mydata[,tail(va.names,1)]
  
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
    Y       <- db$y
    
    
    # output for this loop
    # G, 
    # Gob (the observed part of G)
    # X containing also number of missing links
    # ni is number of students in all schools
    for (i in 1:nsch) {
      cat("School :", i, "/", nsch, "\n")
      schi         <- school[i]
      dbi          <- db[db$sschlcde == schi,]
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
        idxm       <- which(dbi$aid %in% dbi[j, mf_coln])
        idxf       <- which(dbi$aid %in% dbi[j, ff_coln])
        matchim[j] <- length(idxm) # observed male links
        matchif[j] <- length(idxf) # observed female links
        idx        <- c(idxm, idxf)
        Gi[j, idx] <- 1
        
        # missing links
        idxm.miss  <- which(dbi[j, mf_coln] %in% c(mislistmis, -3))
        idxf.miss  <- which(dbi[j, ff_coln] %in% c(mislistmis, -3))
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
      va.zx       <- mydata[c(tmp[x] + 1):tmp[x+1],z]   
      matrix(kronecker(va.zx, va.zx, FUN = dist1), sch.size[x])}))
  }))
  # for va.log2
  X2tmp         <- do.call("cbind", lapply(va.log2, function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx     <- mydata[c(tmp[x] + 1):tmp[x+1],z]   
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

################################ Estimation #############################
######## dataset and model settings
dataset    <- data.frame(GPA = Y, X)
Xlogit     <- as.data.frame(Xlogit)
Model      <- GPA ~ female + hispanic + raceblack + raceasian + raceother + melhigh + memhigh + memiss + mjprof + mjother + mjmiss + age
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

# SGMM estimator
set.seed(123)
obs.sgmm.est     <- smmSAR(formula    = Model,
                           contextual = TRUE,
                           dnetwork   = lapply(G, function(x) 1*(x > 0)), # there is no missing link 
                           smm.ctr    = list(R = 1L, iv.power = 2L, opt.tol = 1e-7, print = TRUE),
                           data       = dataset)
saveRDS(obs.sgmm.est, file = "obs.sgmm.est.RDS")
summary(obs.sgmm.est)

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

# SGMM estimator
# We first estimate the distribution of the network using a logit regression on the observed network
Aobs.mis          <- as.logical(mat.to.vec(lapply(G, function(x) 1*(x > 0))))[as.logical(mat.to.vec(Gobmis))]
Xlogit.mis        <- as.matrix(Xlogit[as.logical(mat.to.vec(Gobmis)),])
miss.logit        <- glm(Aobs.mis ~ Xlogit.mis, family = binomial(link = "logit"), weights = weights)
summary(miss.logit)
rho.mis           <- miss.logit$coefficients
var.rho           <- summary(miss.logit)$cov.unscaled

# estimated distribution for unobserved links
dnetwork.miss     <- vec.to.mat(plogis(c(cbind(1, as.matrix(Xlogit))%*%rho.mis)), N = sch.size)
dnetwork.miss     <- lapply(1:nsch, function(x) (G[[x]] > 0)*Gobmis[[x]] + dnetwork.miss[[x]]*(1 - Gobmis[[x]]))
  
set.seed(123)
miss.sgmm.est     <- smmSAR(formula    = Model,
                            contextual = TRUE,
                            dnetwork   = dnetwork.miss, # there is no missing link 
                            smm.ctr    = list(R = 1000L, iv.power = 2L, opt.tol = 1e-4, print = TRUE),
                            data       = dataset)
saveRDS(miss.sgmm.est, file = "miss.sgmm.est.RDS")
summary(miss.sgmm.est)

# variance computation
fdist.miss        <- function(rho.mis, var.rho, Xlogit, G, Gobmis){
  M               <- length(G)
  N               <- sapply(G, nrow)
  rho.sim         <- MASS::mvrnorm(mu = rho.mis, Sigma = var.rho)
  dsim            <- vec.to.mat(plogis(c(cbind(1, as.matrix(Xlogit))%*%rho.sim)), N = N)
  lapply(1:M, function(x) (G[[x]] > 0)*Gobmis[[x]] + dsim[[x]]*(1 - Gobmis[[x]]))
}
fdist_args        <- list(rho.mis = rho.mis, var.rho = var.rho, Xlogit = Xlogit, G = G, Gobmis = Gobmis)

smiss.sgmm.est     <- summary(miss.sgmm.est, dnetwork = dnetwork.miss, data = dataset, 
                             .fun = fdist.miss, .args = fdist_args, sim = 500L, ncores = 4L)
saveRDS(smiss.sgmm.est, file = "smiss.sgmm.est.RDS")

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
  tcm0  <- matrix(rep((nmx - est.nfriendm[idx])/(nmx - friendm[[x]]), each = sch.size[x]), sch.size[x], byrow = TRUE)
  tcf1  <- matrix(rep(est.nfriendf[idx]/friendf[[x]], each = sch.size[x]), sch.size[x], byrow = TRUE)
  tcf0  <- matrix(rep((nfx - est.nfriendf[idx])/(nfx - friendf[[x]]), each = sch.size[x]), sch.size[x], byrow = TRUE)
  tcm1[is.na(tcm1)] <- 0
  tcm0[is.na(tcm0)] <- 0
  tcf1[is.na(tcf1)] <- 0
  tcf0[is.na(tcf0)] <- 0
  tcm   <- tcm1*G[[x]] + tcm0*(1 - G[[x]])
  tcf   <- tcf1*G[[x]] + tcf0*(1 - G[[x]])
  tcf*wf*gendf[[x]] + tcm*wm*(1 - gendf[[x]])
})
weights <- mat.to.vec(weights)[as.logical(mat.to.vec(Gobtmis))]

# Network distribution
Aobs.tmis         <- as.logical(mat.to.vec(lapply(G, function(x) 1*(x > 0))))[as.logical(mat.to.vec(Gobtmis))]
Xlogit.tmis       <- as.matrix(Xlogit[as.logical(mat.to.vec(Gobtmis)),])
tmiss.logit       <- glm(Aobs.tmis ~ Xlogit.tmis, family = binomial(link = "logit"), weights = weights)
summary(tmiss.logit)
rho.tmis          <- tmiss.logit$coefficients
var.rho           <- summary(tmiss.logit)$cov.unscaled


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

# SGMM estimator
dnetwork.tmiss    <- vec.to.mat(plogis(c(cbind(1, as.matrix(Xlogit))%*%rho.tmis)), N = sch.size)
dnetwork.tmiss    <- lapply(1:nsch, function(x) (G[[x]] > 0)*Gobtmis[[x]] + dnetwork.tmiss[[x]]*(1 - Gobtmis[[x]]))

set.seed(123)
tmiss.sgmm.est    <- smmSAR(formula    = Model,
                            contextual = TRUE,
                            dnetwork   = dnetwork.tmiss, # there is no missing link 
                            smm.ctr    = list(R = 1000L, iv.power = 2L, opt.tol = 1e-4, print = TRUE),
                            data       = dataset)
saveRDS(tmiss.sgmm.est, file = "tmiss.sgmm.est.RDS")
summary(tmiss.sgmm.est)

# variance computation
fdist.tmiss       <- function(rho.tmis, var.rho, Xlogit, G, Gobtmis){
  M               <- length(G)
  N               <- sapply(G, nrow)
  rho.sim         <- MASS::mvrnorm(mu = rho.tmis, Sigma = var.rho)
  dsim            <- vec.to.mat(plogis(c(cbind(1, as.matrix(Xlogit))%*%rho.sim)), N = N)
  lapply(1:M, function(x) (G[[x]] > 0)*Gobtmis[[x]] + dsim[[x]]*(1 - Gobtmis[[x]]))
}
fdist_args        <- list(rho.tmis = rho.tmis, var.rho = var.rho, Xlogit = Xlogit, G = G, Gobtmis = Gobtmis)

stmiss.sgmm.est   <- summary(tmiss.sgmm.est, dnetwork = dnetwork.tmiss, data = dataset, 
                              .fun = fdist.tmiss, .args = fdist_args, sim = 500L, ncores = 4L)
saveRDS(stmiss.sgmm.est, file = "stmiss.sgmm.est.RDS")
print(stmiss.sgmm.est)


####### Plot simulations
obs.ntw           <- obs.bayes.est$posterior
recons.ntw        <- miss.bayes.est$posterior$theta
form.ntw          <- miss.bayes.est$posterior$rho
Xnames            <- c("Female", "Hispanic", paste("Race =", c("Black", "Asian", "Other")),
                       paste0("Mother Edu ", c("< High", "> High", "= Missing")), 
                       paste("Mother Job =", c("Professional", "Other", "Missing")), "Age")

XnamesRho         <- c("Intercept", "Same sex", "Both Hispanic", "Both White", "Both Black",
                       "Both Asian", "Mums Educ < high", "Mums Educ > high",
                       "Mums Job Professional", "Age absolute diff")

# simulations
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

# densities
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

# Network formation model
c1 = "blue"
par(mfrow = c(2, 5), mar = c(2, 2, 2, 2.1))
for (i in 1:10) {
  plot(form.ntw[,i], type = "l", main = XnamesRho[i], col = c1, xlab = "", ylab = "")
}

##################### Policy simulation
# Centrality
library(doParallel)
# observed network
cent.obs    <- unlist(lapply(norm.network(G), function(x) solve(diag(nrow(x)) - colMeans(tail(obs.est$posterior, 1e4))["Gy"]*x, rep(1, nrow(x)))))

# reconstructed network
f.cent.rec  <- function(p.log, alpha){
  network   <- 1*(as.matrix(cbind(1, Xlogit[,c("same.sex", va.log1[-1], "diff.age")])) %*% p.log + rlogis(nrow(Xlogit)) > 0)
  network.n <- vec.to.mat(network, sch.size)
  network.n <- lapply(1:nsch, function(x) Gobs[[x]]*G[[x]] + (1 - Gobs[[x]])*network.n[[x]])
  network.n <- norm.network(network.n)
  unlist(lapply(network.n, function(x) solve(diag(nrow(x)) - alpha*x, rep(1, nrow(x)))))
}
cent.tmiss  <- mclapply(19900:20000, function(x) f.cent.rec(miss.est$posterior$rho[x,], miss.est$posterior$theta[x, "Gy"]), mc.cores = 3)
cent.tmiss  <- apply(as.data.frame(cent.tmiss), 1, mean) 

centdata    <- data.frame(observed = cent.obs, reconstructed = cent.tmiss)

save(cent.tmiss, file = "~/Dropbox/Papers - In progress/Partial Network/AddHealth/centrality.rda")
ggplot(data = centdata, aes(x = observed, y = reconstructed)) + geom_point()

# # Simulate the impact of shock on the network
# # estimate logit on the observed network
nFobs       <- unlist(lapply(G, function(x) rowSums(x > 0)))
nFmis       <- unlist(nmatch)
nFr         <- nFobs + nFmis
mean(nFobs); mean(nFr)
sd(nFobs); sd(nFr)
Xlogit$link <- mat.to.vec(G)

logit0      <- glm(paste0("link ~ ", paste0(c("same.sex", va.log1[-1], "diff.age"), collapse = "+")), family = "binomial", data = Xlogit)
summary(logit0)
# 
# This function simulate impact of policy on the network and y
conterfact  <- function(p.log, p.sar, impact){
  print(impact)
  p.log[1]  <- p.log[1] + impact
  network   <- 1*(as.matrix(cbind(1, Xlogit[,c("same.sex", va.log1[-1], "diff.age")])) %*% p.log + rlogis(nrow(Xlogit)) > 0)
  network.n <- vec.to.mat(network, sch.size, normalise = TRUE)
  ysim      <- CDatanet::simsar(formula    = ~ female + hispanic + raceblack + raceasian + raceother + melhigh + memhigh + memiss + mjprof + mjother + mjmiss + age,
                                contextual = TRUE,
                                Glist      = network.n,
                                theta      = p.sar,
                                data       = mydata)$y
  nfriends  <- unlist(lapply(network.n, function(x) rowSums(x > 0)))
  out       <- c(mean(nfriends), sd(nfriends), mean(ysim), sd(ysim))
  names(out)<- c("mean.nfriends", "sd.nfriends", "mean.y", "sd.y")
  out
}

# application using the observed network
impact      <- c(sort(seq(-3, -0, 0.05), decreasing = TRUE), seq(0.05, 0.5, 0.05))
p.sar.obs   <- apply(tail(obs.est$posterior, 1e4), 2, mean)
p.sar.obs   <- c(p.sar.obs["Gy"], p.sar.obs[names(p.sar.obs) != "Gy"])
p.log.obs   <- logit0$coefficients
policy.obs  <- sapply(impact, function(x) conterfact(p.log.obs, p.sar.obs, x))

# application using the reconstructed network
p.sar.nobs  <- apply(tail(miss.est$posterior$theta, 1e4), 2, mean)
p.sar.nobs  <- c(p.sar.nobs["Gy"], p.sar.nobs[names(p.sar.nobs) != "Gy"])
p.log.nobs  <- apply(tail(miss.est$posterior$rho, 1e4), 2, mean)
policy.nobs <- sapply(impact, function(x) conterfact(p.log.nobs, p.sar.nobs, x))


policy      <- cbind(impact = impact, t(policy.obs), t(policy.nobs))
write.csv(policy, row.names = FALSE, file = "~/Dropbox/Papers - In progress/Partial Network/AddHealth/policy.csv")
