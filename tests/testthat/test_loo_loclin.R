context("Test loo_loclin()")
library(looperr)

n1 = sample(100, 1)
k1 = 3
X1 = matrix(c(rep(1,n1), rnorm(k1*n1, mean=5)), nrow=n1)
Y1 = X1%*%1:(k1+1)
H1 = diag(rep(1, k1))

test_that("loo_loclin Correct output dimensions", {
  expect_equal()
})
