set.seed(1234)
bench = list()
bench_df = list()
results = data.frame()
for(n in c(10, 100, 1000, 10000, 100000, 1000000)){
  for(grps in c(10, 100, 1000, 10000)){
    if(n*grps>10000000){
      next
    }
    # Simulate values
    X = cbind(rep(1, n*grps), rnorm(n*grps))
    w = rep(1, grps*n)
    g = rep(1:grps, each=n)
    y = rnorm(1)*X[,2] + rnorm(n*grps, sd=0.5)

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
      unit="s",
      times=10
    )
    print(bench[[paste0(n, "_", grps)]])
    bench_df[[paste0(n, "_", grps)]] = cbind(
      summary(bench[[paste0(n, "_", grps)]]), list(nObs=n, nGroups=grps)
    )
    results = rbind(results, bench_df[[paste0(n, "_", grps)]])
  }
}
results = results[, c("expr","median", "nObs", "nGroups")]
colnames(results) = c("expr", "seconds", "nObs", "nGroups")
fwrite(results, file="data-raw/results_R.csv")
