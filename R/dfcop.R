#' Discrete factor copula models
#' 
#' Fit discrete factor copula models to data.
#'
#' @param data a binary matrix or data.frame.
#' @param family_set a character vector of families; as in `dfcop_dist()`,
#'   see *Details* for additional options.
#' @param selcrit criterion for family selection, either `"loglik"`, `"aic"`, or 
#'   `"bic"`.
#' @param ngrid number of nodes and weights for the Gaussian quadrature.
#' @param keep_data whether the data should be stored (necessary for computing
#'   fit statistics and using `fitted()`).
#' 
#' @details
#' In addition to the families in `dfcop_dist()`, the following convenience
#' definition from  can be used (and combined):
#' \describe{
#' \item{"all"}{all families}.
#' \item{"archimedean"}{archimedean families.}
#' \item{"elliptical"}{elliptical families.}
#' \item{"bbs"}{BB families.}
#' \item{"oneparametric "}{one parameter families.}
#' \item{"twoparametric "}{two parameter families.}
#' }
#' Partial matching is activated. For example, you can write `"arch"` instead 
#' of the full name `"archimedean"`.
#'
#'
#' @return An object inherting from `dfcop` and `dfcop_dist`.
#'
#' @examples
#' # Sample size, dimension, correlation parameter
#' n <- 1e3
#' d <- 1e1
#' rho <- 0.7
#' 
#' # The model
#' prob <- runif(d)
#' dfcop <- dfcop_dist(prob, "gaussian", rho)
#' 
#' # Data
#' X <- rdfcop(n, dfcop)
#' 
#' # Fit
#' fit <- dfcop(X)
#' fit
#' 
#' @importFrom stats optim
#' @export
dfcop <- function(data, 
                  family_set = "all", 
                  selcrit = "bic",
                  ngrid = 100,
                  keep_data = TRUE) {
  
  # check the data
  stopifnot(is.matrix(data) || is.data.frame(data))
  if(is.data.frame(data)) 
    data <- as.matrix(data)
  stopifnot(ncol(data) >= 2)
  stopifnot(all(data == 0 || data == 1))
  
  # family_set can only use standard family names in cpp
  family_set <- family_set_all_defs[pmatch(family_set, family_set_all_defs)]
  family_set <- expand_family_set(family_set)
  check_family_set(family_set)
  
  # other checks
  selcrit <- match.arg(selcrit, c("aic", "bic", "loglik"))
  stopifnot(length(keep_data) == 1 && (keep_data == TRUE || keep_data == FALSE))
  stopifnot(length(ngrid) == 1 && ngrid > 0)
  
  # fit margins
  prob <- apply(data, 2, mean)
  
  # fit all copula models
  rule <- gauss.quad.prob(ngrid, dist = "uniform", l = 0, u = 1)
  fits <- lapply(family_set, function(family) fitOne(family, prob, data, rule))
  
  # select best model according to selcrit
  lls <- -sapply(fits, function(x) x$nll)
  npars <- sapply(fits, function(x) x$npar)
  if (selcrit == "loglik") {
    best <- which.max(lls)
  } else if (selcrit == "aic") {
    best <- which.min(-2*lls+2*npars)
  } else {
    best <- which.min(-2*lls+log(nrow(data))*npars)
  }
  fit <- dfcop_dist(prob, family_set[best], fits[[best]]$parameters)
  fit$npar <- fits[[best]]$npar
  fit$vcov <- fits[[best]]$vcov
  
  # add the data
  if (keep_data) {
    fit$data <- data
  }
  
  structure(fit, class = c("dfcop", "dfcop_dist"))
}

fitOne <- function(family, prob, data, rule) {
  nll <- function(x, family) {
    pv <- sapply(prob, function(p) 
      hbicop(cbind(p, rule$nodes), cond_var = 2, family, 0, x))
    -sum(log(dfcop_pdf_cpp(rule$weights, pv, data)))
  }
  if (family == "indep") {
    out <- list(nll = -sum(log(ddfcop(data, prob))),
                parameters = numeric(0),
                npar = 0, 
                vcov = 0)
  } else {
    bounds <- if_vec_to_matrix(get_bounds(family))
    npar <- length(bounds[,2])
    tmp <- optim_better(bounds[,2], 
                        function(par) nll(par, family),
                        lower = bounds[,1]+1e-2, upper = bounds[,3]-1e-2, 
                        method = ifelse(npar == 1, "Brent", "L-BFGS-B"),
                        hessian = TRUE)
    if (is.null(tmp$err)) {
      H <- tmp[[1]]$hessian/nrow(data)
      vcov <- ifelse(npar == 1, 1/as.numeric(H), solve(H))
      out <- list(nll = tmp[[1]]$value,
                  parameters = tmp[[1]]$par,
                  npar = npar,
                  vcov = vcov)
    } else {
      out <- list(nll = NA,
                  parameters = NA,
                  npar = NA,
                  vcov = NA)
    }
  }
  return(out)
}

#' S3 methods for discrete factor copula model
#'
#' Currently available:
#' \describe{
#' \item{`predict`}{Predict the probabilities of new observations.}
#' \item{`fitted`}{Fitted probabilities.}
#' \item{`logLik`}{Log-likelihood at the MLE.}
#' \item{`AIC`}{Akaike's Information Criterion (see \link{AIC}).}
#' \item{`BIC`}{Bayesian Information Criterion (see \link{BIC}).}
#' \item{`coef`}{The MLE.}
#' \item{`vcov`}{Variance-covariance matrix at the MLE.}
#' }
#'
#' @aliases fitted.dfcop logLik.dfcop AIC.dfcop BIC.dfcop coef.dfcop vcov.dfcop
#'
#' @param object a `dfcop` object.
#' @param newdata points where the fit shall be evaluated.
#' @param ... unused.
#'
#' @examples
#' # Sample size, dimension, correlation parameter
#' n <- 1e3
#' d <- 1e1
#' rho <- 0.7
#' 
#' # The model
#' prob <- runif(d)
#' dfcop <- dfcop_dist(prob, "gaussian", rho)
#' 
#' # Data
#' X <- rdfcop(n, dfcop)
#' 
#' # Fit
#' fit <- dfcop(X)
#' all.equal(predict(fit, X), fitted(fit))
#' 
#' # Other methods
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#' 
#' @export
predict.dfcop <- function(object, newdata, ...) {
  if (is.null(newdata) && is.null(object$data))
    stop("no newdata and data have not been stored, 
         provide newdata or use keep_data = TRUE when fitting.")
  newdata <- if_vec_to_matrix(newdata)
  ddfcop(object$data, object)
}

#' @rdname predict.dfcop
#' @export
fitted.dfcop <- function(object, ...) {
  if (is.null(object$data))
    stop("data have not been stored, use keep_data = TRUE when fitting.")
  ddfcop(object$data, object)
}

#' @rdname predict.dfcop
#' @export
logLik.dfcop <- function(object, ...) {
  if (is.null(object$data))
    stop("data have not been stored, use keep_data = TRUE when fitting.")
  val <- sum(log(ddfcop(object$data, object)))
  attr(val, "nobs") <- nrow(object$data)
  attr(val, "df") <- object$npar
  class(val) <- "logLik"
  val
}

#' @rdname predict.dfcop
#' @export
coef.dfcop <- function(object, ...) as.numeric(object$bicop$parameters)

#' @rdname predict.dfcop
#' @export
vcov.dfcop <- function(object, ...) object$vcov
