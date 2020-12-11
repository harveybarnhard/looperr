#' @useDynLib looperr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @exportPattern "^[[:alpha:]]+"

.onUnload <- function (libpath) { library.dynam.unload("looperr", libpath)}
