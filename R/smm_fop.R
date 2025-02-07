# The functions starting by fgmm compute the objective function of the gmm estimator
# The functions starting by fbeta compute beta, gamma and se as function of alpha
# the following number 0, 1, 2, 3  stands for Gy and GX are observed, GX is observed and not Gy, 
# Gy is observed and not GX, neither Gy nor GX is observed
# the following nc means no contextual effects
# pr means that the process will be printed


# GX observed, Gy not observed
optim1       <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, V, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  fgmm1(alpha, R, S, T, distr, Ilist, y, X1, X2, V, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum)
}

optim1pr     <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, V, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  out        <- fgmm1(alpha, R, S, T, distr, Ilist, y, X1, X2, V, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum)
  cat("alpha:", alpha, "** objective:", out, "\n")
  out
}

optim1nc     <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  fgmm1nc(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum)
}

optim1ncpr   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  out        <- fgmm1nc(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum)
  cat("alpha:", alpha, "** objective:", out, "\n")
  out
}

optim1fe     <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, V, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  fgmm1fe(alpha, R, S, T, distr, Ilist, y, X1, X2, V, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum)
}

optim1fepr   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, V, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  out        <- fgmm1fe(alpha, R, S, T, distr, Ilist, y, X1, X2, V, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum)
  cat("alpha:", alpha, "** objective:", out, "\n")
  out
}

optim1ncfe   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  fgmm1ncfe(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum)
}

optim1ncfepr <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  out        <- fgmm1ncfe(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum)
  cat("alpha:", alpha, "** objective:", out, "\n")
  out
}


# GX not observed, Gy not observed
optim3   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  fgmm3(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum)
}

optim3pr   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  # print(runif(2))
  out      <- fgmm3(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum)
  cat("alpha:", alpha, "** objective:", out, "\n")
  out
}

optim3nc   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  fgmm3nc(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum)
}

optim3ncpr   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  out        <- fgmm3nc(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum)
  cat("alpha:", alpha, "** objective:", out, "\n")
  out
}

optim3fe     <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  fgmm3fe(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum)
}

optim3fepr   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  out        <- fgmm3fe(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, Kx2, ninstr, M, N, Pm, Ncum)
  cat("alpha:", alpha, "** objective:", out, "\n")
  out
}

optim3ncfe   <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  fgmm3ncfe(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum)
}

optim3ncfepr <- function(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum, seed){
  assign(".Random.seed", seed, envir = .GlobalEnv)#This is used in the function smmSAR. Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function smmSAR
  out        <- fgmm3ncfe(alpha, R, S, T, distr, Ilist, y, X1, X2, W, smoother, hN, Kx1, ninstr, M, N, Pm, Ncum)
  cat("alpha:", alpha, "** objective:", out, "\n")
  out
}