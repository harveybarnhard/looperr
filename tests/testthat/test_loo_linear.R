context("Test loo_loclin()")
library(looperr)

test_that("fastols compared to lm, no weights", {
  set.seed(1234)
  n = 200
  X = cbind(rep(1, n),seq(1,10, length.out=n))
  w = rep(1, n)
  y = sin(X[,2]) + rnorm(n, sd=0.5)
  H = X%*%solve(crossprod(X, X))%*%t(X)
  fastout = fastols(X, y, w)
  origout = lm(y ~ X[,2])
  expect_equal(as.vector(fastout$beta),
               as.vector(origout$coefficients))
  expect_equal(as.vector(fastout$hatdiag),
               as.vector(diag(H)))
  expect_equal(as.vector(fastout$loo_pred_err),
               as.vector(origout$residuals/(1-diag(H))))
})

test_that("fastols compared to lm, with weights", {
  set.seed(1234)
  n = 200
  X = cbind(rep(1, n),seq(1,10, length.out=n))
  w = runif(n, 0, 1)
  W = diag(w)
  y = sin(X[,2]) + rnorm(n, sd=0.5)
  H = X%*%solve(crossprod(X, W)%*%X)%*%t(X)%*%W
  fastout = fastols(X, y, w)
  origout = lm(y ~ X[,2], weights=w)
  expect_equal(as.vector(fastout$beta),
               as.vector(origout$coefficients))
  expect_equal(as.vector(fastout$hatdiag),
               as.vector(diag(H)))
  expect_equal(as.vector(X%*%fastout$beta),
               as.vector(origout$fitted.values))
  expect_equal(as.vector(y-X%*%fastout$beta),
               as.vector(y-origout$fitted.values))
  expect_equal(as.vector(fastout$loo_pred_err),
               as.vector(origout$residuals/(1-diag(H))))
})

