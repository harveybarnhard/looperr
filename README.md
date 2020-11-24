# looperr
An R Package for evaluating leave-one-out (loo) prediction error (perr)
of linear smoothers.
Currently only works for linear regression and local polynomial regression
using a multivariate Gaussian kernel.

Fast implementation in C++ using stable QR decomposition method.

![](examples/looperr_example1.png)

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

One can go further and easily analyze the leave-one-out prediction error point-by-point.
This is often more useful than analyzing residuals because a residual can be small
simply due to that point exerting high-leverage on the line of best fit. But
if you were to leave out that point from the regression, you might actually find
a large residual from the prediction of the response on the left out point.
Here's an example of how to analyze the leave-one-out prediction errors as outputted
from the calls to `linsmooth()` above.

```r
# Smooth the prediction errors
predvals_meta = list()
predvals_meta[[3]] = linsmooth(X, predvals[[3]]$loo_pred_err)
predvals_meta[[4]] = linsmooth(X, predvals[[4]]$loo_pred_err)

# Plot the prediction error points and the fitted curves
plot(x=X[,2], y=predvals[[4]]$loo_pred_err, col="#CC0000", pch=16, cex=0.9)
points(x=X[,2], y=predvals[[3]]$loo_pred_err, col="#00CCCC", pch=16, cex=0.9)
lines(x=X[,2], y=predvals_meta[[4]]$fitted.values, col="#CC0000", lwd=2)
lines(x=X[,2], y=predvals_meta[[3]]$fitted.values, col="#00CCCC", lwd=2)

# Add a legend
legend(x=1,y=-1.4,
       legend=c("Linear Regression" ,
                "Local Linear Regression"),
       col=c("#CC0000", "#00CCCC"), lty=1, lwd=2, box.lty=0, bg=NA)
```

TODO: compare residuals vs. prediction error.

![](examples/looperr_example2.png)

# Setup and Installation

```r
devtools::install_github("harveybarnhard/looperr")
```
