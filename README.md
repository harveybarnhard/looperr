# looperr
<!-- badges: start -->
  [![R build status](https://github.com/harveybarnhard/looperr/workflows/R-CMD-check/badge.svg)](https://github.com/harveybarnhard/looperr/actions)
  [![CRAN Status](https://www.r-pkg.org/badges/version/looperr)](https://www.r-pkg.org/badges/version/looperr)
<!-- badges: end -->

An R package for quickly performing **l**inear sm**oo**thing by grou**p** (**loop**) while
simultaneously computing prediction error (**err**). Or maybe an R package
for quickly performing leave-one-out (**loo**) prediction error (**perr**).

Fast implementation in C++ using OpenMP and non-redundant matrix decomposition methods.

# Smoothing Methods
<table>
    <thead>
        <tr>
            <th>Method</th>
            <th>Submethod</th>
            <th>Multivariate</th>
            <th>By Group?</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td rowspan=2>Least-Squares Regression</td>
            <td rowspan=1>Unweighted</td>
            <td>✔️</td>
            <td>✔️</td>
        </tr>
        <tr>
            <td>Weighted</td>
            <td>✔️</td>
            <td>✔️</td>
        </tr>
        <tr>
            <td rowspan=3>Local Linear Regression</td>
            <td>Uniform Kernel</td>
            <td>❌</td>
            <td>✔️</td>
        </tr>
        <tr>
            <td>Epanechnikov Kernel</td>
            <td>❌</td>
            <td>❌</td>
        </tr>
        <tr>
            <td>Gaussian Kernel</td>
            <td>✔️</td>
            <td>❌</td>
        </tr>
    </tbody>
</table>

# Setup and Installation for Windows, Linux, and Mac

<table>
    <thead>
        <tr>
            <th>OS</th>
            <th>Version Tested</th>
            <th>Build Passing</th>
            <th>Out-of-Package OpenMP Support</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td rowspan=1>Windows</td>
            <td rowspan=1>Microsoft Windows 10.0</td>
            <td>✔️</td>
            <td>✔️ </td>
        </tr>
        <tr>
            <td rowspan=1>Mac</td>
            <td rowspan=1>Mac OS X 10.15.7 (Catalina)</td>
            <td>✔️</td>
            <td>❌ (see extra setup installations)</td>
        </tr>
        <tr>
            <td rowspan=1>Linux</td>
            <td rowspan=1>Ubuntu 20.04</td>
            <td>✔️</td>
            <td>✔️</td>
        </tr>
    </tbody>
</table>

Installation is the same for Windows, Linux, and Mac. Just type the following in the R Console:
```r
devtools::install_github("harveybarnhard/looperr")
```
However, since the default compiler of Mac OS does not support OpenMP,
you will not benefit from the increased efficiency of parallelization
unless you perform some extra setup steps. See the 
[looperr wiki page](https://github.com/harveybarnhard/looperr/wiki/Extra-Setup-and-Installation-for-Mac)
for these extra installation steps.

# Usage
The main function of this package is `linsmooth()` which
performs a similar role to the standard `lm()` function for
regressions. By default, `linsmooth()` performs linear regression.
If your dataset is `X`, your response vector is `y`, and your
group identifier is `g`. Here is a smörgåsbord of
commands you might want to run:

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
# Local linear regression with Gaussian kernel and bandwidth H=1
linsmooth(X, y, method="loclin", H=1)
# Local linear regression with Uniform Kernel, bandwidth H=1, by group g
linsmooth(X, y, method="loclin", H=1, bygroup=g)
```

![](examples/looperr_example1.png)
# Theory
This package uses a bunch of linear algebra "tricks" to reduce
the number of computations, leading to quicker runtimes
and better numerical stability.
