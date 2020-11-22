set.seed(1234)
n = 200
X = cbind(rep(1, n),seq(1,10, length.out=n))
y = sin(X[,2]) + rnorm(n, sd=0.5)
H = list(diag(c(5,5)), diag(c(1,1)), diag(c(0.05,0.05)))
predvals = list()
for(i in 1:3){
  predvals[[i]] = loclin_gauss(X, H[[i]], y)$pred_vals
}
png(filename="looperr_example.png", width=800, height=480)
plot(x=X[,2], y=y,
     axes=FALSE, xaxt="n", yaxt="n", ann=FALSE,
     pch=16, cex=0.9)
lines(X[,2], predvals[[1]], col="#00CCCC", lwd=2)
lines(X[,2], predvals[[2]], col="#FFCC33", lwd=2)
lines(X[,2], predvals[[3]], col="#CC0000", lwd=2)
legend(x=0.5,y=-1,
       legend=c("Big Bandwidth", "Medium Bandwidth", "Small Bandwidth"),
       col=c("#00CCCC", "#FFCC33", "#CC0000"), lty=1, lwd=2, box.lty=0, bg=NA)
dev.off()

