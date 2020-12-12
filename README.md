# looperr
<!-- badges: start -->
  [![R build status](https://github.com/harveybarnhard/looperr/workflows/R-CMD-check/badge.svg)](https://github.com/harveybarnhard/looperr/actions)
  [![CRAN Status](https://www.r-pkg.org/badges/version/looperr)](https://www.r-pkg.org/badges/version/looperr)
<!-- badges: end -->

An R package for quickly performing **l**inear sm**oo**thing by grou**p** (**loop**) while
simultaneously computing prediction error (**err**). Or maybe an R package
for quickly performing leave-one-out (**loo**) prediction error (**perr**).

Smoothing methods include:
* Linear regression
* Weighted linear regression
* local polynomial regression with uniform, Gaussian, or Epanechnikov kernel

Fast implementation in C++ using OpenMP and stable QR methods. 

![](examples/looperr_example1.png)

# Setup and Installation
The package has been tested and runs on Windows 10 and Ubuntu 20.04.
Changes are being made so that the package will also run on Mac OS X.
```r
devtools::install_github("harveybarnhard/looperr")
```

# Usage
The main function of this package is `linsmooth()` which
performs a similar role to the standard `lm()` function for
regressions. The default is to run a local linear regression using
a Gaussian kernel. If your dataset is `X` and your response vector is `y`,
then the following few lines 

```r
# Linear regression
linsmooth(X,y)
# Weighted linear regression
linsmooth(X,y, w=wts)
# Linear regression by group
linsmooth(X,y, g=bygroup)
# Weighted linear regression by group
linsmooth(X,y, w=wts, g=bygroup)
# Local linear regression with Gaussian kernel and optimal bandwidth
linsmooth(X, y, method="loclin")
# Local linear regression with Epanechnikov kernel and optimal bandwidth
linsmooth(X, y, method="loclin", kernel="epan")
# Local linear regression with Epanechnikov kernel and bandwidth H=1
linsmooth(X, y, method="loclin", H=1)
```

# Theory
This package uses a bunch of linear algebra "tricks" to reduce
the number of computations, leading to quicker runtimes
and better numerical stability. For proofs and explanations,
see my
[blog post](https://harveybarnhard.com/posts/evaluating-prediction-error.html)
on this package. The blog post is currently incomplete, and
more proofs will be added later.
