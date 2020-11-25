context("Test loo_loclin()")
library(looperr)

test_that("loo_loclin Correct LOO error", {
  n = 200
  X = cbind(rep(1, n),seq(1,10, length.out=n))
  y = sin(X[,2]) + rnorm(n, sd=0.5)
  h = runif(1, 0.1, 1)
  looind = sample(1:n, size=1)
  total = loclin(X, matrix(1), y, X, 1, kernel=1L)
  loout = loclin(X[-looind,], matrix(1), y[-looind], X, 0, kernel=1L)
  total_looer = total$loo_pred_err[looind]
  loout_looer = y[looind]-loout$fitted.values[looind]
  expect_equal(loout_looer, total_looer)
})

test_that("loo_loclin Correct LOO error", {
  n = 200
  X = cbind(rep(1, n),seq(1,10, length.out=n))
  y = sin(X[,2]) + rnorm(n, sd=0.5)
  h = runif(1, 0.1, 1)
  looind = sample(1:n, size=1)
  total = loclin(X, matrix(1), y, X, 1, kernel=2L)
  loout = loclin(X[-looind,], matrix(1), y[-looind], X, 0, kernel=2L)
  total_looer = total$loo_pred_err[looind]
  loout_looer = y[looind]-loout$fitted.values[looind]
  expect_equal(loout_looer, total_looer)
})

