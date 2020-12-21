context("Test loo_loclin()")
library(looperr)

test_that("fastols compared to lm, no weights", {
  n = 200
  X = cbind(rep(1, n), rnorm(n))
  w = rep(1, n)
  y = rnorm(1)*X[,2] + rnorm(n, sd=0.5)
  H = X%*%solve(crossprod(X, X))%*%t(X)
  fastout = fastols(X, y, compute_se=1, compute_hat=1)
  origout = lm(y ~ X[,2])
  origoutsum = summary(origout)
  expect_equal(as.vector(fastout$beta),
               as.vector(origout$coefficients))
  expect_equal(as.vector(fastout$hatdiag),
               as.vector(diag(H)))
  expect_equal(as.vector(fastout$fitted.values),
               as.vector(origout$fitted.values))
  expect_equal(as.vector(y-fastout$fitted.values),
               as.vector(y-origout$fitted.values))
  expect_equal(as.vector(fastout$loo_pred_err),
               as.vector(origout$residuals/(1-diag(H))))
  expect_equal(fastout$VCV,
               origoutsum$sigma^2*unname(as.matrix(origoutsum$cov.unscaled)))
})

test_that("fastols compared to lm, with weights", {
  n = 200
  X = cbind(rep(1, n), rnorm(n))
  w = runif(n, 0, 1)
  W = diag(w)
  y = rnorm(1)*X[,2] + rnorm(n, sd=0.5)
  H = X%*%solve(crossprod(X, W)%*%X)%*%t(X)%*%W
  fastout = fastolswt(X, y, w, compute_se=1, compute_hat=1)
  origout = lm(y ~ X[,2], weights=w)
  origoutsum = summary(origout)
  expect_equal(as.vector(fastout$beta),
               as.vector(origout$coefficients))
  expect_equal(as.vector(fastout$hatdiag),
               as.vector(diag(H)))
  expect_equal(as.vector(fastout$fitted.values),
               as.vector(origout$fitted.values))
  expect_equal(as.vector(y-fastout$fitted.values),
               as.vector(y-origout$fitted.values))
  expect_equal(as.vector(fastout$loo_pred_err),
               as.vector((y-origout$fitted.values)/(1-diag(H))))
  expect_equal(fastout$VCV,
               origoutsum$sigma^2*unname(as.matrix(origoutsum$cov.unscaled)))
})

test_that("fastols_by compared to lm", {
  n = 1000
  X = cbind(rep(1, n), rnorm(n))
  g = rep(c(1,2,3,4,5), each=200)
  y = rnorm(1)*X[,2] + rnorm(n, sd=0.5)
  fastout = fastols_by(X, y, g, 1, compute_se=1, compute_hat=1)
  for(i in 1:5){
    Xj = X[g==i,]
    H = Xj%*%solve(crossprod(Xj, Xj))%*%t(Xj)
    origout = lm(y[g==i] ~ X[g==i,2])
    origoutsum = summary(origout)
    expect_equal(as.vector(fastout$beta[,i]),
                 as.vector(origout$coefficients))
    expect_equal(as.vector(fastout$hatdiag[g==i]),
                 as.vector(diag(H)))
    expect_equal(y[g==i] - as.vector(fastout$residuals)[g==i],
                 as.vector(origout$fitted.values))
    expect_equal(as.vector(fastout$residuals)[g==i],
                 as.vector(y[g==i]-origout$fitted.values))
    expect_equal(as.vector(fastout$loo_pred_err[g==i]),
                 as.vector((y[g==i]-origout$fitted.values)/(1-diag(H))))
    expect_equal(fastout$variance[,i],
                 origoutsum$sigma^2*unname(diag(origoutsum$cov.unscaled)))
  }
})

test_that("fastols_by compared to lm: two cores", {
  n = 1000
  grps = 1000
  X = cbind(rep(1, n), rnorm(n*grps))
  g = rep(1:grps, each=n)
  y = rnorm(1)*X[,2] + rnorm(n*grps, sd=0.5)
  fastout = fastols_by(X, y, g, 2, compute_se=1, compute_hat=1)
  for(i in sample(1:1000, 10)){
    Xj = X[g==i,]
    H = Xj%*%solve(crossprod(Xj, Xj))%*%t(Xj)
    origout = lm(y[g==i] ~ X[g==i,2])
    origoutsum = summary(origout)
    expect_equal(as.vector(fastout$beta[,i]),
                 as.vector(origout$coefficients))
    expect_equal(as.vector(fastout$hatdiag[g==i]),
                 as.vector(diag(H)))
    expect_equal(y[g==i] - as.vector(fastout$residuals)[g==i],
                 as.vector(origout$fitted.values))
    expect_equal(as.vector(fastout$residuals)[g==i],
                 as.vector(y[g==i]-origout$fitted.values))
    expect_equal(as.vector(fastout$loo_pred_err[g==i]),
                 as.vector((y[g==i]-origout$fitted.values)/(1-diag(H))))
    expect_equal(fastout$variance[,i],
                 origoutsum$sigma^2*unname(diag(origoutsum$cov.unscaled)))
  }
})

test_that("fastols_by compared to lm: two cores + more vars", {
  n = 1000
  grps = 1000
  X = cbind(rep(1, n), rnorm(n*grps), rnorm(n*grps))
  g = rep(1:grps, each=n)
  y = rnorm(1)*X[,2] + rnorm(1)*X[,3] + rnorm(n*grps, sd=0.5)
  fastout = fastols_by(X, y, g, 2, compute_se=1, compute_hat=1)
  for(i in sample(1:1000, 10)){
    Xj = X[g==i,]
    H = Xj%*%solve(crossprod(Xj, Xj))%*%t(Xj)
    origout = lm(y[g==i] ~ X[g==i,2:3])
    origoutsum = summary(origout)
    expect_equal(as.vector(fastout$beta[,i]),
                 as.vector(origout$coefficients))
    expect_equal(as.vector(fastout$hatdiag[g==i]),
                 as.vector(diag(H)))
    expect_equal(y[g==i] - as.vector(fastout$residuals)[g==i],
                 as.vector(origout$fitted.values))
    expect_equal(as.vector(fastout$residuals)[g==i],
                 as.vector(y[g==i]-origout$fitted.values))
    expect_equal(as.vector(fastout$loo_pred_err[g==i]),
                 as.vector((y[g==i]-origout$fitted.values)/(1-diag(H))))
    expect_equal(fastout$variance[,i],
                 origoutsum$sigma^2*unname(diag(origoutsum$cov.unscaled)))
  }
})

test_that("fastols_by compared to lm: two cores + weights", {
  n = 1000
  grps = 1000
  X = cbind(rep(1, n), rnorm(n*grps))
  w = runif(n*grps, 0, 1)
  g = rep(1:grps, each=n)
  y = rnorm(1)*X[,2] + rnorm(n*grps, sd=0.5)
  fastout = fastols_bywt(X, y, w, g, 2, compute_se=1, compute_hat=1)
  for(i in sample(1:1000, 10)){
    Xj = X[g==i,]
    Wj = diag(w[g==i])
    H = Xj%*%solve(crossprod(Xj, Wj)%*%Xj)%*%t(Xj)%*%Wj
    origout = lm(y[g==i] ~ X[g==i,2], weights=w[g==i])
    origoutsum = summary(origout)
    expect_equal(as.vector(fastout$beta[,i]),
                 as.vector(origout$coefficients))
    expect_equal(as.vector(fastout$hatdiag[g==i]),
                 as.vector(diag(H)))
    expect_equal(y[g==i] - as.vector(fastout$residuals)[g==i],
                 as.vector(origout$fitted.values))
    expect_equal(as.vector(fastout$residuals)[g==i],
                 as.vector(y[g==i]-origout$fitted.values))
    expect_equal(as.vector(fastout$loo_pred_err[g==i]),
                 as.vector((y[g==i]-origout$fitted.values)/(1-diag(H))))
    expect_equal(fastout$variance[,i],
                 origoutsum$sigma^2*unname(diag(origoutsum$cov.unscaled)))
  }
})

