context("Test loo_loclin()")
library(looperr)

set.seed(1234)
n = 200
X = cbind(rep(1, n),seq(1,10, length.out=n))
y = sin(X[,2]) + rnorm(n, sd=0.5)

test_that("loo_loclin Correct output dimensions", {
  expect_equal()
})

