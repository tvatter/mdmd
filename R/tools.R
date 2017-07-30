#' Internal: Turn vector input into a matrix with two columns
#'
#' @param u input data
#'
#' @return either a matrix with two columns, or an error if u is neither a
#' matrix, data.frame, or a length two vector
#'
#' @noRd
if_vec_to_matrix <- function(u) {
  if (NCOL(u) == 1)
    u <- matrix(u, 1, length(u))
  u
}

#' Internal: Convert arguments to `dfcop_dist` object.
#' @param prob the marginal probabilities as passed in function call.
#' @param family the family as passed in function call.
#' @param parameters the parameters as passed in function call.
#' @return A `dfcop_dist` object.
#' @noRd
args2dfcop <- function(prob, family = "indep", parameters = numeric(0)) {
  if (all(inherits(prob, "dfcop_dist"))) {
    return(prob)
  } else {
    return(dfcop_dist(prob, family, parameters))
  }
}

factory <- function(fun)
  function(...) {
    warn <- err <- NULL
    res <- withCallingHandlers(
      tryCatch(fun(...), error=function(e) {
        err <<- conditionMessage(e)
        NULL
      }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(res, warn=warn, err=err)
  }

optim_better <- factory(optim)