#' @title The PartialNetwork package
#' @description The \pkg{PartialNetwork} package implements instrumental variables (IV) and Bayesian estimators for the linear-in-mean SAR model (e.g. Bramoulle et al., 2009) when
#' the distribution of the network is available, but not the network itself. To make the computations faster \pkg{PartialNetwork} uses \code{C++} through the \pkg{Rcpp} package (Eddelbuettel et al., 2011). 
#'
#' @details 
#' Two main functions are provided to estimate the linear-in-mean SAR model using only the distribution of the network. The function
#' \code{\link{sim.IV}} generates valid instruments using the distribution of the network (see Propositions 1 and 2 in Boucher and Houndetoungan (2025)). Once the instruments are constructed,
#' one can estimate the model using standard IV estimators. We recommend the function \link[AER]{ivreg} 
#' from the package \pkg{AER} (Kleiber et al., 2020). The function \link{mcmcSAR} performs a Bayesian estimation based on an adaptive MCMC (Atchade and Rosenthal, 2005). In that case, 
#' the distribution of the network acts as prior distribution for the network.\cr
#' The package \pkg{PartialNetwork} also implements a network formation model based on Aggregate Relational Data (McCormick and Zheng, 2015; Breza et al., 2017). This part of the package
#' relies on the functions \link{rvMF}, \link{dvMF} and \link{logCpvMF} partly implemented in C++, but using code from \pkg{movMF} (Hornik and Grun, 2014).
#' 
#' @references 
#' Atchade, Y. F., & Rosenthal, J. S., 2005, On adaptive markov chain monte carlo algorithms, \emph{Bernoulli}, 11(5), 815-828, \doi{10.3150/bj/1130077595}.
#' @references 
#' Boucher, V., & Houndetoungan, A., 2025. Estimating peer effects using partial network data. \emph{arXiv preprint arXiv:2509.08145}.
#' @references 
#' Bramoulle, Y., Djebbari, H., & Fortin, B., 2009, Identification of peer effects through social networks, \emph{Journal of econometrics}, 150(1), 41-55, \doi{10.1016/j.jeconom.2008.12.021}.
#' @references Breza, E., Chandrasekhar, A. G., McCormick, T. H., & Pan, M., 2020, Using aggregated relational data to feasibly identify network structure without network data, \emph{American Economic Review}, 110(8), 2454-84, \doi{10.1257/aer.20170861}
#' @references Eddelbuettel, D., Francois, R., Allaire, J., Ushey, K., Kou, Q., Russel, N., ... & Bates, D., 2011,
#' \pkg{Rcpp}: Seamless \R and \code{C++} integration, \emph{Journal of Statistical Software}, 40(8), 1-18, \doi{10.18637/jss.v040.i08}
#' @references  
#' Lee, L. F., 2004, Asymptotic distributions of quasi-maximum likelihood estimators for spatial autoregressive models. Econometrica, 72(6), 1899-1925, \doi{10.1111/j.1468-0262.2004.00558.x}
#' @references 
#' LeSage, J. P. 1997, Bayesian estimation of spatial autoregressive models, \emph{International regional science review}, 20(1-2), 113-129, \doi{10.1177/016001769702000107}.
#' @references 
#' Mardia, K. V., 2014, Statistics of directional data, \emph{Academic press}.
#' @references McCormick, T. H., & Zheng, T., 2015, Latent surface models for networks using Aggregated Relational Data, 
#' \emph{Journal of the American Statistical Association}, 110(512), 1684-1695, \doi{10.1080/01621459.2014.991395}.
#' @references  
#' Wood, A. T., 1994, Simulation of the von Mises Fisher distribution. \emph{Communications in statistics-simulation and computation}, 23(1), 157-164. \doi{10.1080/03610919408813161}.
#' @useDynLib PartialNetwork, .registration = TRUE
"_PACKAGE"
NULL
