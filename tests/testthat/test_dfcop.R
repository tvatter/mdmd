context("Class 'dfcop'")

test_that("dummy call creates proper dfcop object", {
  X <- matrix(rbinom(1e3*3,1,0.5), ncol = 3)
  fit <- dfcop(X, "indep")
  expect_s3_class(fit, "dfcop_dist")
  expect_s3_class(fit, "dfcop")
  expect_identical(names(fit), c("prob", "bicop", "npar", "vcov", "data"))
})

test_that("family sets works", {
  expect_error(check_family_set("imnotafamily"))
  sapply(family_set_defs, expand_family)
})

test_that("aic/bic for dfcop works", {
  set.seed(0)
  # Sample size, dimension, correlation parameter
  n <- 1e3
  d <- 1e1
  rho <- 0.7
  
  # The model
  prob <- runif(d)
  dfcop <- dfcop_dist(prob, "gaussian", rho)
  
  # Data
  X <- rdfcop(n, dfcop)

  # AIC/BIC select the proper family
  expect_identical(dfcop(X, selcrit = "aic")$bicop$family, "gaussian")
  expect_identical(dfcop(X, selcrit = "bic")$bicop$family, "gaussian")

})

test_that("S3 methods for dfcop work", {
  set.seed(0)
  # Sample size, dimension, correlation parameter
  n <- 1e3
  d <- 1e1
  rho <- 0.7
  
  # The model
  prob <- runif(d)
  dfcop <- dfcop_dist(prob, "gaussian", rho)
  
  # Data
  X <- rdfcop(n, dfcop)
  
  # Fit
  fit <- dfcop(X)
  
  # The methods work as expected
  expect_equal(predict(fit, X), fitted(fit))
  expect_s3_class(logLik(fit), "logLik")
  expect_equivalent(names(attributes(logLik(fit))), c("nobs", "df", "class"))
  expect_is(vcov(fit), "numeric")
  expect_is(coef(fit), "numeric")
  
  # They also throw the right errors
  expect_error(predict(fit))
  fit$data <- NULL
  expect_error(fitted(fit))
  expect_error(logLik(fit))
})