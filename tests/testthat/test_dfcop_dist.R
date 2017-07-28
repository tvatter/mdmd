context("Class 'dfcop_dist'")

test_that("constructor creates proper dfcop_dist object", {
  dist <- dfcop_dist(c(0.5, 0.5), "gumbel", 3)
  expect_s3_class(dist, "dfcop_dist")
  expect_identical(names(dist), c("prob", "bicop"))
})

test_that("d/r functions work", {
  set.seed(123)
  
  # For d function
  n <- 1e4
  prob <- runif(10)
  dfcop <- dfcop_dist(prob, "gaussian", 0.5)
  x <- sapply(prob, function(p) rbinom(n, 1, p))
  p <- ddfcop(x, dfcop)
  p2 <- apply(sapply(seq_along(prob), function(i) 
    dbinom(x[,i], 1, prob[i])), 1, prod)
  
  # proba in [0,1]
  expect_gte(min(p), 0)
  expect_lte(max(p), 1)
  
  # integrate to 1
  expect_lte(abs(mean(p/p2) - 1), 1e-2)
  
  # For r function 
  prob <- runif(2)
  dfcop <- dfcop_dist(prob, "gaussian", 0.5)
  x <- rdfcop(n, dfcop)
  x2 <- expand.grid(c(0,1), c(0,1))
  p <- apply(x2, 1, function(y) 
    mean(x[,1] == y[1] & x[,2] == y[2]) - ddfcop(y, dfcop))
  
  # Data in {0,1}
  expect_true(all(x == 0 || x == 1))
  
  # Empirical probabilities close to theoretical ones
  expect_true(all(abs(p) < 1e-1))
})