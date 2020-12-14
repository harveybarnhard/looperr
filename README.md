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
            <td rowspan=1>Microsoft Windows Server 2019 10.0.17763</td>
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
## Extra Setup and Installation for Mac (Speed Improvements)
If you have a Mac and want to take advantage of parallelization
for fast run-times, then you'll likely have to take some extra steps
if your default compiler (clang on Mac) does not have OpenMP support.
If you're on a Windows or Linux machine, no need to read further--you're
already benefitting from OpenMP parallelization.

The easiest way to get OpenMP support for looperr on a Mac 
is to change the default R compiler to gcc from
clang. Here's a step-by-ste on how to do this via the Mac Terminal.
First, we will need to install OpenMP.
```shell
brew install libomp
```
If you have any problems, you might have to install the 
[Homebrew package manager](https://brew.sh/).
Next, we want to install our new OpenMP enabled compiler gcc. Before installing,
we remove gfortran (fortran compiler) because `brew install gcc` will install
its own version
```shell
rm '/usr/local/bin/gfortran'
brew install gcc
```
Finally, we create a Makevars file which sets the default C, C++, and Fortran
compilers (in this case to gcc version 10). You can create this file "by hand"
but you might as well just create it in the terminal as follows:
```shell
mkdir ~/.R
touch ~/.R/Makevars
{
  echo VER=-10
  echo CC=gcc-10
  echo CXX=g++-10
  echo CXX11=g++-10
  echo CXX14=g++-10
  echo CXX17=g++-10
  echo FC=/usr/bin/gfortran
  echo F77=/usr/bin/gfortran
  echo CFLAGS=-mtune=native -g -O2 -Wall -pedantic
  echo CXXFLAGS=-mtune=native -g -O2 -Wall -pedantic
  echo F77="gfortran-4.8"
  echo FC="gfortran-4.8"
  echo FLIBS = ""
} >~/.R/Makevars
```

# Usage
The main function of this package is `linsmooth()` which
performs a similar role to the standard `lm()` function for
regressions. By default, `linsmooth()` performs linear regression.
If your dataset is `X`, your response vector is `y`, and your
group identifier is `g` the here is a smörgåsbord of
command commands you might want to run:

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
and better numerical stability. For proofs and explanations,
see my
[blog post](https://harveybarnhard.com/posts/evaluating-prediction-error.html)
on this package. The blog post is currently incomplete, and
more proofs will be added later.
