# looperr
An R Package for evaluating leave-one-out (loo) prediction error (perr)
of linear smoothers.
Currently only works for linear regression and local polynomial regression
using a multivariate Gaussian kernel.

Fast implementation in C++ using stable QR decomposition method.

![](examples/looperr_example.png)

# Theory

# Example

The figure at the top of this README was created using the following code:

```r
# Create some fake sinusoidal data
set.seed(1234)
n = 200
X = cbind(rep(1, n),seq(1,10, length.out=n))
y = sin(X[,2]) + rnorm(n, sd=0.5)

# Perform local linear regression with three bandwidth matrices
predvals = list()
predvals[[1]] = linsmooth(X, y, H=matrix(2))    # Large bandwidth
predvals[[2]] = linsmooth(X, y)                 # No bandwidth => choose optimal
predvals[[3]] = linsmooth(X, y, H=matrix(0.1))  # Small bandwidth

# Plot the local linear fits using base R graphics
plot(x=X[,2], y=y,
     axes=FALSE, xaxt="n", yaxt="n", ann=FALSE,
     pch=16, cex=0.9)
lines(X[,2], predvals[[1]]$fitted.values, col="#00CCCC", lwd=2)
lines(X[,2], predvals[[2]]$fitted.values, col="#FFCC33", lwd=2)
lines(X[,2], predvals[[3]]$fitted.values, col="#CC0000", lwd=2)

# Add a legend
legend(x=0.5,y=-1,
       legend=c("Large Bandwidth: 2" ,
                paste0("Optimal Bandwidth: ",
                       round(predvals[[2]]$bandwidth,2)),
                "Small Bandwidth: 0.1"),
       col=c("#00CCCC", "#FFCC33", "#CC0000"), lty=1, lwd=2, box.lty=0, bg=NA)
```
