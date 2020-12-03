set.seed(1234)
bench = list()
results = data.frame()
for(n in c(1000, 10000, 100000, 1000000)){
  for(grps in c(5, 50, 100, 200, 500)){
    # Simulate values
    X = cbind(rep(1, n), rnorm(n))
    w = rep(1, n)
    g = rep(1:grps, each=n/grps)
    y = rnorm(1)*X[,2] + rnorm(n, sd=0.5)

    # Create dataframe and save
    outdf = data.frame(constant = X[,1], x=X[,2], w=w, g=g, y=y)
    fwrite(outdf, file=paste0("data-raw/benchmark", n,"_", grps,".csv"))

    # Benchmark in R
    message("nObs=", n," nGroups=", grps)
    bench[[paste0(n, "_", grps)]] = microbenchmark::microbenchmark(
      baseR = for(i in 1:grps){
        lm(y[g==i] ~ X[g==i,2])
      },
      cpp_onecore = fastols_by(X,y,w,g,1),
      cpp_twocore = fastols_by(X,y,w,g,2),
      cpp_threecore = fastols_by(X,y,w,g,3),
    unit="s")
    print(bench[[paste0(n, "_", grps)]])
  }
  bench[[paste0(n, "_", grps)]] = cbind(
    summary(bench[[paste0(n, "_", grps)]]), list(nObs=n, nGroups=grps)
  )
  results = rbind(results, summary(bench[[1]]))
}
