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

# Smooth the residuals =========================================================
resid_diff = y-predvals[[3]]$fitted.values - y-predvals[[4]]$fitted.values
resid_smooth = linsmooth(X, resid_diff)
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

# 3D example ===================================================================
dt = data.table::data.table(volcano)
dt[, Y:=1:nrow(volcano)]
dt = data.table::melt(dt, measure.vars=patterns("^V"), value.name="Z",
                      variable.factor=FALSE)
data.table::setnames(dt, old="variable", new="X")
dt[, X:= as.numeric(substr(X, 2, nchar(X)))]
dt[, Znoise := Z + rnorm(nrow(dt), mean=0, sd=10)]
test = linsmooth(cbind(dt$X, dt$Y), dt$Znoise, H=c(1,2))
dt[, Zsmooth := linsmooth(X, )]

fig = plot_ly(dt, x= ~X, y= ~Y, z= ~Znoise)
fig = fig %>% add_markers(size=1, color=)
add_markers(plot_ly(dt, x= ~X, y= ~Y, z= ~Znoise), size=1)

# Blog post example ============================================================
library(pBrackets)
library(ggplot2)
library(latex2exp)

# Simulate data
set.seed(123)
n = 20
X = rnorm(n, mean=0, sd=0.2)
beta = 0.5
Y = X*beta + rnorm(n, mean=0, sd=0.1)

# Create a high-leverage point
Y[18] = Y[18] + .5
X[18] = X[18] + .2

# Plot the regression
ggplot(data.frame(X, Y), aes(x=X, y=Y)) +
        geom_point() +
        geom_smooth(method="lm", formula="y~x", se=FALSE,
                    aes(color="With High-Leverage Point")) +
        geom_smooth(data=data.frame(X=X[-18], Y=Y[-18]),
                    method="lm", formula="y~x", se=FALSE,
                    aes(x=X, y=Y,color="Without High-Leverage Point"),
                    fullrange=TRUE) +
        coord_cartesian(xlim=c(min(X), max(X))) +
        geom_curve(aes(x = X[18]+0.1, y = Y[18]-0.1,
                       xend = X[18]+0.01, yend = Y[18]-0.008),
                   colour = "#CC0000",
                   size=1,
                   curvature = 0.1,
                   arrow = arrow(length = unit(0.03, "npc"))) +
        annotate(geom="text",x=X[18]+0.11, y=Y[18]-0.11, label="High Leverage Point",
                 size=6) +
        annotate(geom="text",x=-0.225, y=0.137,
                 label=latex2exp::TeX("$\\hat{\\epsilon}_i$"),
                 size=6, parse=TRUE) +
        annotate(geom="text",x=-0.146 , y=0.1067,
                 label=latex2exp::TeX("$\\hat{\\epsilon}_i^{leaveout}$"),
                 size=6, parse=TRUE) +
        scale_colour_manual(name="", values=c("#CC0000", "#00CCCC")) +
        theme_void() +
        theme(legend.position=c(0.8, 0.15),
              legend.text=element_text(size=16))
grid.brackets(158, 474, 158, 32)
grid.brackets(165, 32, 165, 546)
