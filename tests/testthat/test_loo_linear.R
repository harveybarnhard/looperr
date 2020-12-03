context("Test loo_loclin()")
library(looperr)

test_that("fastols compared to lm, no weights", {
  n = 200
  X = cbind(rep(1, n), rnorm(n))
  w = rep(1, n)
  y = rnorm(1)*X[,2] + rnorm(n, sd=0.5)
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
  n = 200
  X = cbind(rep(1, n), rnorm(n))
  w = runif(n, 0, 1)
  W = diag(w)
  y = rnorm(1)*X[,2] + rnorm(n, sd=0.5)
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

test_that("fastols_by compared to lm", {
  n = 1000
  X = cbind(rep(1, n), rnorm(n))
  w = rep(1, n)
  g = rep(c(1,2,3,4,5), each=200)
  y = rnorm(1)*X[,2] + rnorm(n, sd=0.5)
  fastout = fastols_by(X, y, w, g)
  for(i in 1:5){
    expect_equal(as.vector(fastout$beta[,i]),
                 as.vector(coefficients(lm(y[g==i] ~ X[g==i,2]))))
  }
})

