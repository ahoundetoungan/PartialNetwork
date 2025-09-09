############################################################################################################
################ This script imports AddHealth data and prepares them for the estimation ###################
############################################################################################################
rm(list = ls())
library(dplyr)
library(PartialNetwork)
library(haven)

InDataPath  <- "A/B/C/PATH_TO_DATAIN/" # Where Add Health data are saved (/ at the end is important)
OutDataPath <- "A/B/C/PATH_TO_DATAOUT/" # Where prepared data for each outcome are saved (/ at the end is important)

# Importing data sets
# Friendship data set (WAVE I)
sfriend  <- read_xpt(paste0(InDataPath, "sfriend.xpt")) %>% arrange(SQID) %>% 
  filter(!(SQID %in% c("999999", ""))) %>%  # Remove student with missing questionnaire ID (cannot be matched)
  mutate(across(ends_with("AID"), as.character)) 

# Inschool data set (WAVE I)
Inschool <- read_xpt(paste0(InDataPath, "Inschool.xpt")) %>% arrange(SQID) %>% 
  filter(!(SQID %in% c("999999", "")), AID != "") %>% # Remove student with missing questionnaire ID and missing student ID (cannot be matched)
  mutate(across(c("SQID", "AID", "SSCHLCDE"), as.character))

# Variable to keep in the Inschool data set
to_keep  <- c("SQID", "AID", "SSCHLCDE", "S1", "S2", "S4", "S6A", "S6B", "S6C", "S6E", "S10A", "S10B", "S10C", "S10D", 
              "S12", "S14")
Inschool <- Inschool %>% select(all_of(to_keep)) 

# Add the friendship data set to the Inschool data set
Inschool <- Inschool %>% left_join(sfriend, by = "SQID")


# Create variables that will be used
Inschool <- Inschool %>% 
  # I recoded some variables before constructing the ones that will be used in the models.
  # Replace 99 with NA
  mutate(across("S1", ~ ifelse(. == 99, NA, .))) %>%
  # Replace 9 with NA
  mutate(across("S2", ~ ifelse(. == 9, NA, .))) %>%
  # Replace values out of 1:4 with NA (this is for the GPA)
  mutate(across(all_of(c("S10A", "S10B", "S10C", "S10D")), ~ ifelse(. %in% (1:4), ., NA))) %>%
  # Model variable creation
  mutate(age = S1, age2 = age^2, S1 = NULL, # age
         male = as.integer(S2 == 1), female = as.integer(S2 == 2), S2 = NULL, # gender
         hispanic = ifelse(S4 %in% c(8, 9), 0, S4), S4 = NULL, # hispanic, replace do not know and multiple responses with No
         racewhite = as.integer(S6A == 1), S6A = NULL, # race
         raceblack = as.integer(S6B == 1), S6B = NULL,
         raceasian = as.integer(S6C == 1), S6C = NULL, 
         raceother = as.integer(S6E == 1), S6E = NULL,
         gpa = 5 - rowMeans(select(., c("S10A", "S10B", "S10C", "S10D")), na.rm = TRUE), # gpa
         across(all_of(c("S10A", "S10B", "S10C", "S10D")), ~ NULL), 
         # mother education
         melhigh = as.integer(S12 %in% c(1, 2, 9, 10)), # less than HS
         mehigh = as.integer(S12 %in% c(3, 4)), # HS
         memhigh = as.integer(S12 %in% (5:8)), # more than HS
         memiss = as.integer(S12 >= 11 | is.na(S12)),
         S12 = NULL, 
         # mother job
         mjhome = as.integer(S14 %in% c(1, 16:18)), # at home
         mjprof = as.integer(S14 %in% c(2, 3)), # professional
         mjother = as.integer(S14 %in% c(4:15, 19)), # other
         mjmiss = as.integer(S14 >= 20 | is.na(S14))) %>% # missing
  arrange(SSCHLCDE, AID) %>% mutate(SCID = SSCHLCDE)

# error codes
mislist      <- c(55555555, 77777777, 88888888, 99999999, 99959995)

# friend variables
mf_coln      <- paste0("MF", 1:5, "AID")
ff_coln      <- paste0("FF", 1:5, "AID")
f_coln       <- c(mf_coln, ff_coln)

# Check if mislist is ID
if (sum(Inschool$AID %in% mislist) > 0) {
  stop("mislist is an ID")
} else {
  cat("mislist is not an ID: OK", "\n")
}

# list of variable (excluding reference variables for identification)
va.names      <- c("female", "hispanic", "raceblack",   "raceasian",   "raceother", "melhigh", "memhigh",
                   "memiss", "mjprof", "mjother",  "mjmiss", "age", "gpa")

# Are there NAs?
Inschool %>% summarise(across(all_of(va.names), ~ sum(is.na(.))))

# Keep row without NA
Inschool        <- Inschool %>% filter(if_all(all_of(va.names), ~ !is.na(.)))

# remove friend from different groups
# remove self friendship
# remove friend non found
N       <- nrow(Inschool)
dscf    <- rep(0,N)
sfre    <- rep(0,N)
nffr    <- rep(0,N)
for (i in 1:N) {
  for (j in f_coln) {
    k   <- which(Inschool$AID == Inschool[[j]][i])
    # remove if different school
    if (length(k) != 0) { # if friend found
      if(Inschool$SCID[i] != Inschool$SCID[k]) {
        Inschool[[j]][i] <- -1
        dscf[i]          <- dscf[i] + 1
      }
      # remove if self friendship
      if(Inschool$AID[i] == Inschool$AID[k]) {
        Inschool[[j]][i] <- -2
        sfre[i]          <- sfre[i] + 1
      }
    } else { # if friend not found
      if (!((Inschool[[j]][i] %in% mislist) | is.na(Inschool[[j]][i]))) { # if not missing or error code
        Inschool[[j]][i] <- -3
        nffr[i]          <- nffr[i] + 1
      }
    }
  }
}

cat("remove", sum(dscf), "link(s) because students from different schools: their code are recode as -1", "\n")
cat("remove", sum(sfre), "self-friendship(s): their code are recoded as -2", "\n")
cat("remove", sum(nffr), "non-found friends: their code are recoded as -3", "\n")

# Keep schools < School size max
Inschool   <- Inschool %>% group_by(SCID) %>% mutate(SchSize = n()) %>% ungroup() %>%
  filter(SchSize <= 200)
school     <- unique(Inschool$SCID)
nsch       <- length(school)
sch.size   <- (Inschool %>% group_by(SCID) %>% summarize(SchSize = n()))$SchSize
  
# dependent variable
Inschool$y <- Inschool[[tail(va.names,1)]]

mislistmis <- c(55555555, 99999999, 99959995)


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
  Y       <- db[["y"]]
  
  
  # output for this loop
  # G, 
  # Gobs1 and Gobs2 (the observed part of G)
  # X containing also number of missing links
  # ni is number of students in all schools
  for (i in 1:nsch) {
    cat("School :", i, "/", nsch, "\n")
    dbi          <- db[db$SCID == school[i],]
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
      idxm       <- which(dbi$AID %in% unlist(dbi[j, mf_coln]))
      idxf       <- which(dbi$AID %in% unlist(dbi[j, ff_coln]))
      idx        <- c(idxm, idxf)
      matchim[j] <- length(idxm) # observed male links
      matchif[j] <- length(idxf) # observed female links
      Gi[j,idx]  <- 1
      
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
    
    # unmatched
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
tmp     <- gen.data(Inschool)
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
    va.zx     <- Inschool[[z]][c(tmp[x] + 1):tmp[x+1]]   
    matrix(kronecker(va.zx, va.zx, FUN = dist1), sch.size[x])}))
}))
# for va.log2
X2tmp         <- do.call("cbind", lapply(va.log2, function(z) {
  mat.to.vec(lapply(1:nsch, function(x) {
    va.zx     <- Inschool[[z]][c(tmp[x] + 1):tmp[x+1]]   
    matrix(kronecker(va.zx, va.zx, FUN = dist2), sch.size[x])}))
}))
Xlogit        <- cbind(X1tmp, X2tmp)  
colnames(Xlogit) <- c("same.sex", va.log1[-1], "diff.age")

# Objects to be saved
SaObj  <- c("friendf", "friendm", "nmatchf", "nmatchm", "matchf", "matchm", 
            "G", "gendf", "dscf", "Gobmis", "Gobtmis", "Gobtop", "Inschool", 
            "nffr", "nsch", "OutDataPath", "InDataPath", "sch.size", "school", 
            "sfre",  "va.log1", "va.log2", "va.names", "X", "Xlogit", "Y")
save(list = SaObj, file = paste0(OutDataPath, "DataPartialNetwork.rda"))
