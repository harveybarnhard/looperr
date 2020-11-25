# Github example ===============================================================
# Create some fake sinusoidal data
set.seed(1234)
n = 200
X = cbind(rep(1, n),seq(1,10, length.out=n))
y = sin(X[,2]) + rnorm(n, sd=0.5)
H = list(matrix(2), matrix(0.1))

# Perform local linear regression with three bandwidth matrices
predvals = list()
predvals[[1]] = linsmooth(X, y, H=2)    # Large bandwidth
predvals[[2]] = linsmooth(X, y, H=0.1)  # Small bandwidth
predvals[[3]] = linsmooth(X, y)         # Optimal bandwidth

# Perform linear regression
predvals[[4]] = linsmooth(X, y, method="ols")

# Plot the local linear fits using base R graphics
png(filename="examples/looperr_example1.png", width=800, height=480)
plot(x=X[,2], y=y,
     axes=FALSE, xaxt="n", yaxt="n", ann=FALSE,
     pch=16, cex=0.9)
lines(X[,2], predvals[[1]]$fitted.values, col="#FFCC33", lwd=2)
lines(X[,2], predvals[[2]]$fitted.values, col="#CC0000", lwd=2)
lines(X[,2], predvals[[3]]$fitted.values, col="#00CCCC", lwd=2)
legend(x=0.5,y=-1,
       legend=c("Large Bandwidth: 2" ,
                paste0("Optimal Bandwidth: ",
                       round(predvals[[3]]$bandwidth,2)),
                "Small Bandwidth: 0.1"),
       col=c("#FFCC33", "#00CCCC", "#CC0000"), lty=1, lwd=2, box.lty=0, bg=NA)
dev.off()

# Smooth the prediction errors
predvals_meta = list()
predvals_meta[[3]] = linsmooth(X, predvals[[3]]$loo_pred_err)
predvals_meta[[4]] = linsmooth(X, predvals[[4]]$loo_pred_err)

# Plot the prediction error points and the fitted curves
png(filename="examples/looperr_example2.png", width=800, height=480)
plot(x=X[,2], y=predvals[[4]]$loo_pred_err, col="#CC0000", pch=16, cex=0.9,
     ylab = "Leave-One-Out Prediction Error", xlab="")
points(x=X[,2], y=predvals[[3]]$loo_pred_err, col="#00CCCC", pch=16, cex=0.9)
lines(x=X[,2], y=predvals_meta[[4]]$fitted.values, col="#CC0000", lwd=2)
lines(x=X[,2], y=predvals_meta[[3]]$fitted.values, col="#00CCCC", lwd=2)

# Add a legend
legend(x=1,y=-1.4,
       legend=c("Linear Regression" ,
                "Local Linear Regression"),
       col=c("#CC0000", "#00CCCC"), lty=1, lwd=2, box.lty=0, bg=NA)
dev.off()

# Smooth the residuals
resid_loc = y-predvals[[3]]$fitted.values
resid_lin = y-predvals[[4]]$fitted.values
# Plot the prediction error points and the fitted curves
png(filename="examples/looperr_example3.png", width=800, height=480)
plot(x=resid_lin, y=predvals[[4]]$loo_pred_err, col="#CC0000", pch=16, cex=0.9,
     ylab = "Leave-One-Out Prediction Error", xlab="")
points(x=resid_loc, y=predvals[[3]]$loo_pred_err, col="#00CCCC", pch=16, cex=0.9)

# Add a legend
legend(x=1,y=-1.4,
       legend=c("Linear Regression" ,
                "Local Linear Regression"),
       col=c("#CC0000", "#00CCCC"), lty=1, lwd=2, box.lty=0, bg=NA)
dev.off()
