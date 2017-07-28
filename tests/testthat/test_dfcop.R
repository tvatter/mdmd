context("Class 'dfcop'")

test_that("dummy call creates proper dfcop object", {
  X <- matrix(rbinom(1e3*3,1,0.5), ncol = 3)
  fit <- dfcop(X, "indep")
  expect_s3_class(fit, "dfcop_dist")
  expect_s3_class(fit, "dfcop")
  expect_identical(names(fit), c("prob", "bicop", "npar", "data"))
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

  expect_identical(dfcop(X, selcrit = "aic")$bicop$family, "gaussian")
  expect_identical(dfcop(X, selcrit = "bic")$bicop$family, "gaussian")
})