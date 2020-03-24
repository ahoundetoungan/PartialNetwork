#' @title The PartialNetwork package
#' @description The \pkg{PartialNetwork} package implements instrumental variables (IV) and Bayesian estimators for the linear-in-mean SAR model (e.g. Bramoullé et al., 2009) when
#' the distribution of the network is availaible, but not the network itself. To make the computations faster \pkg{PartialNetwork} uses \code{C++} through the \pkg{Rcpp} package (Eddelbuettel et al., 2011). 
#'
#' @details 
#' Two main functions are provided to estimate the linear-in-mean SAR model using only the distribution of the network. The function
#' \code{\link{sim.IV}} generates valid instruments using the distribution of the network (see Propositions 1 and 2 in Boucher and Houndetoungan (2020)). Once the instruments are constructed,
#' one can estimate the model using standard IV estimators. We recommend the function \link[AER]{ivreg} 
#' from the package \pkg{AER} (Kleiber et al., 2020). The function \link{mcmcSAR} performs a Bayesian estimation based on an adaptive MCMC (Atchadé and Rosenthal, 2005). In that case, 
#' the distribution of the network acts as prior distribution for the network.\cr
#' The package \pkg{PartialNetwork} also implements a network formation model based on Aggregate Relational Data (McCormick and Zheng, 2015; Breza et al., 2017). This part of the package
#' relies on the functions \link{rvMF}, \link{dvMF} and \link{logCpvMF} partly implemented in C++, but using code from \pkg{movMF} (Hornik and Grün, 2014).
#' 
#' @references 
#' Atchadé, Y. F., & Rosenthal, J. S. (2005). On adaptive markov chain monte carlo algorithms. \emph{Bernoulli}, 11(5), 815-828. \url{https://projecteuclid.org/euclid.bj/1130077595}.
#' @references 
#' Boucher, V., & Houndetoungan, A. (2020). Estimating peer effects using partial network data. \emph{Draft avaliable at} \url{https://houndetoungan.wixsite.com/aristide/research}.
#' @references 
#' Bramoullé, Y., Djebbari, H., & Fortin, B. (2009). Identification of peer effects through social networks. \emph{Journal of econometrics}, 150(1), 41-55. \url{https://www.sciencedirect.com/science/article/abs/pii/S0304407609000335}.
#' @references Breza, E., Chandrasekhar, A. G., McCormick, T. H., & Pan, M. (2020). Using aggregated relational data to feasibly
#'  identify network structure without network data. \emph{American Economic Review, forthcoming.}, \url{https://arxiv.org/abs/1703.04157}
#' @references Eddelbuettel, D., François, R., Allaire, J., Ushey, K., Kou, Q., Russel, N., ... & Bates, D. (2011),
#' \pkg{Rcpp}: Seamless \R and \code{C++} integration. \emph{Journal of Statistical Software}, 40(8), 1-18.
#' \url{http://www.jstatsoft.org/v40/i08/}.
#' @references Hornik, K., & Grün, B. (2014). \pkg{movMF}: An R package for fitting mixtures of von Mises-Fisher distributions. \emph{Journal of Statistical Software}, 58(10), 1-31.
#' \url{https://epub.wu.ac.at/4893/}.
#' @references 
#' Kleiber, C., Zeileis, A., & Zeileis, M. A. (2020). Package ‘AER’. \R package version 1.2, 4. \url{https://cran.r-project.org/web/packages/AER/index.html}.
#' @references 
#' Mardia, K. V. (2014). Statistics of directional data. Academic press. \url{https://www.elsevier.com/books/statistics-of-directional-data/mardia/978-0-12-471150-1}.
#' @references McCormick, T. H., & Zheng, T. (2015). Latent surface models for networks using Aggregated Relational Data. 
#' Journal of the American Statistical Association, 110(512), 1684-1695. \url{https://www.tandfonline.com/doi/abs/10.1080/01621459.2014.991395}.
#' @references  
#' Wood, A. T. (1994). Simulation of the von Mises Fisher distribution. Communications in statistics-simulation and computation, 23(1), 157-164. \url{https://www.tandfonline.com/doi/abs/10.1080/03610919408813161}.
#' @useDynLib PartialNetwork, .registration = TRUE
"_PACKAGE"
NULL
