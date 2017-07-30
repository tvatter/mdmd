#' A discrete factor copula distribution
#' 
#' A S3 class to store discrete factor copula distributions.
#'
#' @param prob the marginal probabilities, a vector of numbers in (0,1) or a 
#' list of vectors with numbers in (0,1).
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
  if(is.vector(prob)) {
    stopifnot(all(prob > 0 && prob < 1))
    nmax <- rep(2, length(prob))
    binary <- TRUE
  }else if(is.list(prob)){
    lapply(prob, function(x) stopifnot(is.vector(prob) && 
                                         all(x >= 0 && sum(x) != 1)))
    nmax <- sapply(prob, length)
    if(all(nmax==2)){ 
      ## recast into binary dfcop_dist
      prob <- sapply(prob, function(x) x[2])
    }else{
      binary <- FALSE
    }
  }else stop()
  bicop <- bicop_dist(family, 0, parameters)
  dist <- list(prob = prob,
               bicop = bicop,
               binary = binary,
               nmax = nmax)
  structure(dist, class = "dfcop_dist")
}

#' @export
print.dfcop_dist <- function(x, ...) {
  if (x$bicop$family == "indep") {
    cat("Discrete factor copula ('dfcop_dist'): ",
        "dimension = ", length(x$prob),
        ", family = ", x$bicop$family,
        ", binary = ", as.character(x$binary),
        sep = "")
  } else {
    cat("Discrete factor copula ('dfcop_dist'): ",
        "dimension = ", length(x$prob),
        ", family = ", x$bicop$family,
        ", dependence parameters = ", x$bicop$parameters,
        ", binary = ", as.character(x$binary),
        sep = "")
  }
}

#' @rdname dfcop_dist
#' @param x an integer vector of the same length as prob, or a matrix with 
#' number of columns equal to the length of prob
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
  stopifnot((is.vector(x) || is.matrix(x)))
  if(dfcop$binary) stopifnot(all(x == 0 || x == 1))
  else stopifnot(all(sapply(1:length(dfcop$prob), 
                            function(i) x[i] %in% 0:dfcop$nmax[i])))
  x <- if_vec_to_matrix(x)
  stopifnot(ncol(x) == length(dfcop$prob))
  stopifnot(length(ngrid) == 1 && ngrid > 0)
  rule <- gauss.quad.prob(ngrid, dist = "uniform", l = 0, u = 1)
  if(dfcop$binary){
    pv <- sapply(dfcop$prob, function(p) hbicop(cbind(p, rule$nodes), 
                                                cond_var = 2, 
                                                dfcop$bicop))
    pdf <- dfcop_pdf_cpp(rule$weights, pv, x)
    # pvt <- t(pv)
    # apply(x, 1, function(xi) 
    #   sum(apply(pvt^as.vector(xi)*(1-pvt)^(1-as.vector(xi)), 
    #             2, prod)*rule$weights))
  }else{
    pv <- sapply(seq_along(prob), function(i){
      pis <- dfcop$prob[[i]]
      if(x[i]==0){ 
        p2 <- 1
        p1 <- sum(pis[-1])
      }else if(x[i]==dfcop$nmax[i]){
        p2 <- rev(pis)[1]
        p1 <- 0
      }else{
        p2 <- 1 - sum(pis[1:(x[i])])
        p2 <- sum(pis[(1+x[i]):length(pis)])
      }
      hbicop(cbind(p2, rule$nodes), cond_var = 2, dfcop$bicop) - 
        hbicop(cbind(p1, rule$nodes), cond_var = 2, dfcop$bicop)
    })
    sum(apply(pv, 1, prod)*rule$weights)
  }
  return(pdf)
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
  U <- apply(U, 2, function(u) 
    hbicop(cbind(u,V), 2, dfcop$bicop, inverse = TRUE))
  if(dfcop$binary){
    x <- sapply(seq_along(prob), function(i) ifelse(U[,i] >= prob[i], 0, 1))
  }else{
    x <- sapply(seq_along(prob), function(i)
      dfcop$nmax[i] - which.first(U[,i] <= cumsum(rev(prob[[i]]))))
  }
  return(x)
}
