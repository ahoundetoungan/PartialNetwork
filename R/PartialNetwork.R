#' @title The PartialNetwork package
#' @description The \pkg{PartialNetwork} package implements IV Compute IV and Bayesian estimator for linear-in-mean SAR model when
#' only the network distribution is availaible. To make these computations faster \pkg{PartialNetwork} relies on a \code{C++} using \pkg{Rcpp} (Eddelbuettel et al., 2011). 
#'
#' @details 
#' ---- Vincent General detail on the package is here ----
#' @references Breza, E., Chandrasekhar, A. G., McCormick, T. H., & Pan, M. (2017). Using aggregated relational data to feasibly
#'  identify network structure without network data (No. w23491). National Bureau of Economic Research.
#' @references Eddelbuettel, D., François, R., Allaire, J., Ushey, K., Kou, Q., Russel, N., ... & Bates, D. (2011),
#' \pkg{Rcpp}: Seamless \R and \code{C++} integration. \emph{Journal of Statistical Software}, 40(8), 1-18.
#' \url{http://www.jstatsoft.org/v40/i08/}.
#' @references McCormick, T. H., & Zheng, T. (2015). Latent surface models for networks using Aggregated Relational Data. 
#' Journal of the American Statistical Association, 110(512), 1684-1695.
#' @useDynLib PartialNetwork, .registration = TRUE
"_PACKAGE"
NULL
