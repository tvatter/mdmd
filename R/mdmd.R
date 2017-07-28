#' Multivariate Discrete Mixture Distributions
#'
#' @name mdmd
#' @docType package
#' @useDynLib mdmd
#' @importFrom Rcpp evalCpp
#' 
#' @author Thibault Vatter, Damien Ackerer
#
#' @keywords package
#'
#' @examples
#' # A 10-dimensional Gaussian factor copula
#' dfcop <- dfcop_dist(runif(10),"gaussian", 0.5)
#' str(dfcop_dist(runif(10), "gauss", 0.5))
#' 
#' # Evaluate the probability mass function
#' ddfcop(cbind(c(1,1,0,0), c(1,0,1,0)), c(0.5, 0.5), "indep")
#' ddfcop(rbinom(10, 1, 0.5), dfcop)
#' 
#' # Sample data
#' X <- rdfcop(1e3, dfcop)
#' 
#' # Fit
#' fit <- dfcop(X)
NULL
