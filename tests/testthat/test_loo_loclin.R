context("Test loo_loclin()")
library(looperr)

test_that("loo_loclin Correct LOO error, Gaussian Kernel", {
  n = 200
  X = cbind(rep(1, n),seq(1,10, length.out=n))
  y = sin(X[,2]) + rnorm(n, sd=0.5)
  h = runif(1, 0.1, 1)
  looind = sample(1:n, size=1)
  Xloo  = matrix(X[looind,], nrow=1, ncol=2)
  total = loclin_sameX(X, y, matrix(h), kernel=1L)
  loout = loclin_diffX(X[-looind,], y[-looind], matrix(h), Xloo, kernel=1L)
  total_looer = total$loo_pred_err[looind]
  loout_looer = y[looind]-as.vector(loout$fitted.values)
  expect_equal(loout_looer, total_looer)
})

test_that("loo_loclin Correct LOO error, Gaussian Kernel", {
  n = 200
  X = cbind(rep(1, n),seq(1,10, length.out=n))
  y = sin(X[,2]) + rnorm(n, sd=0.5)
  h = runif(1, 0.1, 1)
  looind = sample(1:n, size=1)
  Xloo  = matrix(X[looind,], nrow=1, ncol=2)
  total = loclin_sameX(X, y, matrix(h), kernel=2L)
  loout = loclin_diffX(X[-looind,], y[-looind], matrix(h), Xloo, kernel=2L)
  total_looer = total$loo_pred_err[looind]
  loout_looer = y[looind]-as.vector(loout$fitted.values)
  expect_equal(loout_looer, total_looer)
})

test_that("loo_loclin Correct LOO error, Uniform Kernel", {
  n = 200
  X = cbind(rep(1, n),seq(1,10, length.out=n))
  y = sin(X[,2]) + rnorm(n, sd=0.5)
  h = runif(1, 0.1, 1)
  looind = sample(1:n, size=1)
  Xloo  = matrix(X[looind,], nrow=1, ncol=2)
  total = loclin_sameX(X, y, matrix(h), kernel=3L)
  loout = loclin_diffX(X[-looind,], y[-looind], matrix(h), Xloo, kernel=3L)
  total_looer = total$loo_pred_err[looind]
  loout_looer = y[looind]-as.vector(loout$fitted.values)
  expect_equal(loout_looer, total_looer)
})
