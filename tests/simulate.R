# Simulate values
n=1000; grps=1000
X = cbind(rep(1, n*grps), rnorm(n*grps))
w = rep(1, grps*n)
g = rep(1:grps, each=n)
y = rnorm(1)*X[,2] + rnorm(n*grps, sd=0.5)
ord = order(g, X[,2])
X = X[ord,]
g = g[ord]
y = y[ord]

microbenchmark::microbenchmark(crossprod(X,X))
microbenchmark(fastols_by(X,y,g, nthr=8), loclin_sameX_byunif(X,y,g, 1, nthr=8), times=10)
