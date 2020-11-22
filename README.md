# looperr
An R Package for evaluating leave-one-out (loo) prediction error (perr).
Currently only works for local polynomial regression using a multivariate
Gaussian kernel and an arbitrary positive-definite bandwidth matrix.

Fast implementation in C++ using stable QR decomposition method.

# Use

```r
# Create some fake data
n = 200
X = cbind(rep(1, n),seq(1,10, length.out=n))
y = sin(X[,2]) + rnorm(n, sd=0.5)
H1 = diag(3,3)
H2 = diag(1,1)
H3 = diag(0.1,0.1)

# Perform local linear regression with three bandwidth matrices
sinsmooth1 = loclin_gauss(X, H1, y)
sinsmooth2 = loclin_gauss(X, H2, y)
sinsmooth3 = loclin_gauss(X, H3, y)

# Plot the local linear fits
plot(X[,2], y)
lines(X[,2], sinsmooth1$pred_vals, col="red")
lines(X[,2], sinsmooth2$pred_vals, col="blue")
lines(X[,2], sinsmooth3$pred_vals, col="green")
```
