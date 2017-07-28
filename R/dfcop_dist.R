#' A discrete factor copula distribution
#' 
#' A dicrete factor copula distribution is specified by:
#'
#' @param prob the marginal probabilities, a vector of numbers in (0,1).
#' @param family the copula family, a string containing the family name (see
#' *Details* for all possible families).
#' @param parameters a vector or matrix of copula paramters.
#'
#' @note
#' The evaluation functions can optionally be used with a `dfcop_dist` object,
#' e.g., `ddfcop(c(1, 0), dfcop_dist(c(0.5, 0.5),"indep"))`.
#'
#' @details
#' The implemented families listed below. Partial matching is activated, i.e.,
#' `"gauss"` is equivalent to `"gaussian"`.
#' \describe{
#' \item{`indep`}{Independence copula.}
#' \item{`gaussian`}{Gaussian copula.}
#' \item{`t`}{Student t copula.}
#' \item{`clayton`}{Clayton copula.}
#' \item{`gumbel`}{Gumbel copula.}
#' \item{`frank`}{Frank copula.}
#' \item{`joe`}{Joe copula.}
#' \item{`bb1`}{BB1 copula.}
#' \item{`bb6`}{BB6 copula.}
#' \item{`bb7`}{BB7 copula.}
#' \item{`bb8`}{BB8 copula.}
#' }
#'
#'
#' @return An object of class `dfcop_dist`.
#'
#' @examples
#' # A 10-dimensional Gaussian factor copula
#' dfcop_dist(runif(10),"gaussian", 0.5)
#' str(dfcop_dist(runif(10), "gauss", 0.5))
#'
#' # A 100-dimensional rotated Clayton factor copula
#' dfcop <- dfcop_dist(runif(100), "clayton", 3)
#' @importFrom rvinecopulib bicop_dist
#' @export
dfcop_dist <- function(prob, 
                       family = "indep", 
                       parameters = numeric(0)) {
  stopifnot(is.vector(prob) && all(prob > 0 && prob < 1))
  bicop <- bicop_dist(family, 0, parameters)
  dist <- list(prob = prob,
               bicop = bicop)
  structure(dist, class = "dfcop_dist")
}

#' @export
print.dfcop_dist <- function(x, ...) {
  if (x$bicop$family == "indep") {
    cat("Discrete factor copula ('dfcop_dist'): ",
        "dimension = ", length(x$prob),
        ", family = ", x$bicop$family,
        sep = "")
  } else {
    cat("Discrete factor copula ('dfcop_dist'): ",
        "dimension = ", length(x$prob),
        ", family = ", x$bicop$family,
        ", dependence parameters = ", x$bicop$parameters,
        sep = "")
  }
}

#' @rdname dfcop_dist
#' @param x a binary vector of the same length as prob
#' @param ngrid number of nodes and weights for the Gaussian quadrature
#' @examples  
#' # evaluate the probability mass function
#' ddfcop(cbind(c(1,1,0,0), c(1,0,1,0)), c(0.5, 0.5), "indep")
#' ddfcop(rbinom(100, 1, 0.5), dfcop)
#' 
#' @importFrom rvinecopulib hbicop
#' @importFrom statmod gauss.quad.prob
#' @export
ddfcop <- function(x, prob, 
                   family = "indep", 
                   parameters = numeric(0),
                   ngrid = 100) {
  dfcop <- args2dfcop(prob, family, parameters)
  stopifnot((is.vector(x) || is.matrix(x)) && all(x == 0 || x == 1))
  x <- if_vec_to_matrix(x)
  stopifnot(ncol(x) == length(dfcop$prob))
  stopifnot(length(ngrid) == 1 && ngrid > 0)
  rule <- gauss.quad.prob(ngrid, dist = "uniform", l = 0, u = 1)
  pv <- sapply(dfcop$prob, function(p) hbicop(cbind(p, rule$nodes), 
                                              cond_var = 2, 
                                              dfcop$bicop))
  dfcop_pdf_cpp(rule$weights, pv, x)
  # pvt <- t(pv)
  # apply(x, 1, function(xi) 
  #   sum(apply(pvt^as.vector(xi)*(1-pvt)^(1-as.vector(xi)), 
  #             2, prod)*rule$weights))
}


#' @rdname dfcop_dist
#' @param n number of observations. 
#' @examples
#' # simulate data
#' rdfcop(1e2, c(0.5, 0.5))
#' rdfcop(1e2, dfcop)
#' @importFrom stats runif
#' @importFrom rvinecopulib rbicop hbicop
#' @export
rdfcop <- function(n, prob, 
                   family = "indep", 
                   parameters = numeric(0)) {
  stopifnot(length(n) == 1 && is.numeric(n))
  n <- round(n)
  dfcop <- args2dfcop(prob, family, parameters)
  prob <- dfcop$prob
  V <- runif(n)
  U <- matrix(runif(n*length(prob)), nrow = n)
  U <- apply(U, 2, function(u) hbicop(cbind(u,V), 2, dfcop$bicop, inverse = TRUE))
  sapply(seq_along(prob), function(i) ifelse(U[,i] >= prob[i], 0, 1))
}
