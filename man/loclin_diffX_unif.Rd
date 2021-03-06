% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{loclin_diffX_unif}
\alias{loclin_diffX_unif}
\title{Function that performs local linear regression
using Gaussian or Epanechnikov kernel.}
\usage{
loclin_diffX_unif(X, y, H, Xeval)
}
\arguments{
\item{X}{an nxk data matrix}

\item{y}{The nx1 output vector}

\item{H}{a kxk positive definite bandwidth matrix}

\item{Xeval}{an mxp matrix at which to predict using local linear regression}
}
\description{
Function that performs local linear regression
using Gaussian or Epanechnikov kernel.
}
