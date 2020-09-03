library(readstata13)
library(ggplot2)
library(PartialNetwork)
library(Matrix)
library(DataCombine)
library(stargazer)
set.seed(123)


# Data
mydata       <- read.dta13("~/tmpah/SchoolData/cleandta.dta")
mydata       <- mydata[order(mydata$sschlcde),]

mislist      <- c(55555555, 77777777, 88888888, 99999999, 99959995)
f_coln       <- c(paste0("mf", 1:5, "aid"), paste0("ff", 1:5, "aid"))


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

# Update data
mydata       <- mydata[mydata$sschlcde %in% school,]
nsch         <- length(school)

# dependent variable
mydata$y     <- mydata[,tail(va.names,1)]

mislistmis   <- c(55555555, 99999999, 99959995)


# This function prepares the data and the network 
gen.data   <- function(db) {
  G       <- vector("list", nsch)
  P       <- vector("list", nsch)
  nmatch  <- vector("list", nsch)
  Y       <- db$y
  
  
  # output for this loop
  # G, 
  # X containing also number of missing links
  # ni is number of students in all schools
  for (i in 1:nsch) {
    cat("School :", i, "/", nsch, "\n")
    schi         <- school[i]
    dbi          <- db[db$sschlcde == schi,]
    Ni           <- nrow(dbi)
    Gi           <- matrix(0, Ni, Ni)
    
    nmatchi      <- numeric() #will contain missing links
    for (j in 1:Ni) {
      idx        <- which(dbi$aid %in% dbi[j, f_coln])
      Gi[j,idx]  <- 1
      
      # missing links
      idx.miss   <- which(dbi[j, f_coln] %in% c(mislistmis, -3))
      nmatchi[j] <- length(idx.miss) # Number of missing links
    }
    
    #  store G
    diag(Gi)     <- 0
    G[[i]]       <- Gi
    
    # store X
    X            <- as.matrix(db[, va.names[-length(va.names)]])    
    
    # Store Y
    Y            <- db$y
    
    # unmatched
    nmatch[[i]]  <- nmatchi
    
    
    # P, the prior distribution of the network
    Pi           <- Gi
    for (j in 1:Ni) {
      if (nmatchi[j] > 0) {
        nfriend              <- sum((Gi[j,] > 0))     
        proba                <- nmatchi[j]/(Ni - 1 - nfriend)
        Pi[j, (Gi[j,] == 0)] <- proba
      }
    }
    diag(Pi) <- 0
    P[[i]]   <- Pi
  }
  
  list(G       = G,
       P       = P,
       X       = X,
       nmatch  = nmatch,
       Y       = Y)
  
}


# Use the function to prepare the data
tmp     <- gen.data(mydata)
G       <- tmp$G
P       <- tmp$P
X       <- tmp$X
nmatch  <- tmp$nmatch
Y       <- tmp$Y


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
n.friend  <- unlist(lapply(G, rowSums))
sum(n.friend)

# number of unmatched friends
nmatchall <- do.call("c", nmatch)
sum(nmatchall)

# number of friends (observed and unmatched)
sum(n.friend) + sum(nmatchall)

# graph missing links
# This reproduces the Figure 4
ggplot(data = data.frame(mm = nmatchall), aes(x = mm)) +
  geom_bar(color = "black", fill = "#eeeeee") + 
  theme_bw() + xlab("") + ylab("Frequency") + 
  scale_x_discrete(limits = 0:10) 
# size 7 × 4 inch


################################ Estimation #############################
########  Observed network as given
Kv         <- 2*ncol(X) + 1  # Number of exogenous explanatory variables

# Hyperparameters (G is the prior distribution of G)
hyperp1    <- list("dnetwork"      = G,
                   "mutheta"       = rep(0,Kv),
                   "invstheta"     = diag(Kv)/100,
                   "muzeta"        = 0,
                   "invszeta"      = 2,
                   "a"             = 4.2,
                   "b"             = (4.2 - 2)*0.5)

obs.est    <- mcmcSAR(formula = Y ~ X,
                      contextual = TRUE,
                      hyperparms = hyperp1,
                      iteration = 10000)


########  Reconstructed network
hyperp2    <- list("dnetwork"      = P,
                   "mutheta"       = rep(0,Kv),
                   "invstheta"     = diag(Kv)/100,
                   "muzeta"        = 0,
                   "invszeta"      = 2,
                   "a"             = 4.2,
                   "b"             = (4.2 - 2)*0.5)

rec.est    <- mcmcSAR(formula    = Y ~ X,
                      contextual = TRUE,
                      hyperparms = hyperp2,
                      iteration  = 10000,
                      ctrl.mcmc  = list(print.level = 2))


####### Summarize results
summary(obs.est)
summary(rec.est)

####### Plot simulations
obs.ntw           <- obs.est$posterior
recons.ntw        <- rec.est$posterior
Xnames            <- c("Female", "Hispanic", paste("Race =", c("Black", "Asian", "Other")),
                       paste0("Mother Edu ", c("< High", "> High", "= Missing")), 
                       paste("Mother Job =", c("Professional", "Other", "Missing")), "Age")

  
### observed network
par(fig = c(0, 1/6, 1 - 0.92/5, 1), mar = c(2, 2, 2, 2.1))
plot(obs.ntw[,26], type = "l", main = "Peer Effect", col = "blue", xlab = "", ylab = "")

par(fig = c(1/6, 2/6, 1 - 0.92/5, 1), new = TRUE)
plot(obs.ntw[,1], type = "l", main = "Intercept", col = "blue", xlab = "", ylab = "")

par(fig = c(2/6, 3/6, 1 - 0.92/5, 1), new = TRUE)
plot(obs.ntw[,27], type = "l", main = expression(sigma^2), col = "blue", xlab = "", ylab = "")

for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (1 + ceiling(i/6))*0.92/5 - 0.04, 1 - ceiling(i/6)*0.92/5 - 0.04),
      new = TRUE)
  plot(obs.ntw[,i + 1], type = "l", main = Xnames[i], col = "blue", xlab = "", ylab = "")
}

for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (3 + ceiling(i/6))*0.92/5 - 0.07999999, 1 - (2 + ceiling(i/6))*0.92/5 - 0.08),
      new = TRUE)
  plot(obs.ntw[,i + 13], type = "l", main = Xnames[i], col = "blue", xlab = "", ylab = "")
}

par(fig = c(0, 1, 1 - 2*0.92/5 - 0.04, 1 - 0.92/5), new = TRUE)
plot(x = 0, main = "Own Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

par(fig = c(0, 1, 1 - 4*0.92/5 - 0.08, 1 - 3*0.92/5 - 0.04), new = TRUE)
plot(x = 0, main = "Contextual Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')


### Reconstructed network
par(fig = c(0, 1/6, 1 - 0.92/5, 1), mar = c(2, 2, 2, 2.1))
plot(recons.ntw[,26], type = "l", main = "Peer Effect", col = "blue", xlab = "", ylab = "")

par(fig = c(1/6, 2/6, 1 - 0.92/5, 1), new = TRUE)
plot(recons.ntw[,1], type = "l", main = "Intercept", col = "blue", xlab = "", ylab = "")

par(fig = c(2/6, 3/6, 1 - 0.92/5, 1), new = TRUE)
plot(recons.ntw[,27], type = "l", main = expression(sigma^2), col = "blue", xlab = "", ylab = "")

for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (1 + ceiling(i/6))*0.92/5 - 0.04, 1 - ceiling(i/6)*0.92/5 - 0.04),
      new = TRUE)
  plot(recons.ntw[,i + 1], type = "l", main = Xnames[i], col = "blue", xlab = "", ylab = "")
}


for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (3 + ceiling(i/6))*0.92/5 - 0.07999999, 1 - (2 + ceiling(i/6))*0.92/5 - 0.08),
      new = TRUE)
  plot(recons.ntw[,i + 13], type = "l", main = Xnames[i], col = "blue", xlab = "", ylab = "")
}

par(fig = c(0, 1, 1 - 2*0.92/5 - 0.04, 1 - 0.92/5), new = TRUE)
plot(x = 0, main = "Own Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

par(fig = c(0, 1, 1 - 4*0.92/5 - 0.08, 1 - 3*0.92/5 - 0.04), new = TRUE)
plot(x = 0, main = "Contextual Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')


# plot densities
par(fig = c(0, 1/6, 1 - 0.92/5, 1), mar = c(2, 2, 2, 2.1))
tmpo <- density(obs.ntw[2001:10000,26])
tmpr <- density(recons.ntw[2001:10000,26])
plot(tmpo, type = "l", main = "Peer Effect",
     col = "blue", xlab = "", ylab = "", xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))))
lines(tmpr, col = "red")

par(fig = c(1/6, 2/6, 1 - 0.92/5, 1), new = TRUE)
tmpo <- density(obs.ntw[2001:10000,1])
tmpr <- density(recons.ntw[2001:10000,1])
plot(tmpo, type = "l", main = "Intercept", col = "blue", xlab = "", ylab = "", 
     xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))),
     ylim = c(min(c(tmpo$y, tmpr$y)), max(c(tmpo$y, tmpr$y))))
lines(tmpr, col = "red")

par(fig = c(2/6, 3/6, 1 - 0.92/5, 1), new = TRUE)
tmpo <- density(obs.ntw[2001:10000,27])
tmpr <- density(recons.ntw[2001:10000,27])
plot(tmpo, type = "l", main = expression(sigma^2), col = "blue", xlab = "", ylab = "",
     xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))),
     ylim = c(min(c(tmpo$y, tmpr$y)), max(c(tmpo$y, tmpr$y))))
lines(tmpr, col = "red")

for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (1 + ceiling(i/6))*0.92/5 - 0.04, 1 - ceiling(i/6)*0.92/5 - 0.04),
      new = TRUE)
  tmpo <- density(obs.ntw[2001:10000, i + 1])
  tmpr <- density(recons.ntw[2001:10000, i + 1])
  plot(tmpo, type = "l", main = Xnames[i], col = "blue", xlab = "", ylab = "",
       xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))),
       ylim = c(min(c(tmpo$y, tmpr$y)), max(c(tmpo$y, tmpr$y))))
  lines(tmpr, col = "red")
}


for (i in 1: 12) {
  par(fig = c(((i - 1) %% 6)/6, ((i - 1) %% 6 + 1)/6,
              1 - (3 + ceiling(i/6))*0.92/5 - 0.07999999, 1 - (2 + ceiling(i/6))*0.92/5 - 0.08),
      new = TRUE)
  tmpo <- density(obs.ntw[2001:10000, i + 13])
  tmpr <- density(recons.ntw[2001:10000, i + 13])
  plot(tmpo, type = "l", main = Xnames[i], col = "blue", xlab = "", ylab = "",
       xlim = c(min(c(tmpo$x, tmpr$x)), max(c(tmpo$x, tmpr$x))),
       ylim = c(min(c(tmpo$y, tmpr$y)), max(c(tmpo$y, tmpr$y))))
  lines(tmpr, col = "red")
}

par(fig = c(0.5, 1, 0.75, 1), new = TRUE)
plot(x = 0, y = 0, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
legend(0, 1.2, legend=c("Observed Network", "Reconstructed Network"),
       col=c("blue", "red"), lty=1:2, cex=1, box.lty=0)

par(fig = c(0, 1, 1 - 2*0.92/5 - 0.04, 1 - 0.92/5), new = TRUE)
plot(x = 0, main = "Own Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

par(fig = c(0, 1, 1 - 4*0.92/5 - 0.08, 1 - 3*0.92/5 - 0.04), new = TRUE)
plot(x = 0, main = "Contextual Effects", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
