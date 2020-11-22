#' This function performs linear smoothing
#'
#' @param X numeric matrix; an nxk data matrix
#' @param y numeric vector;an nx1 response vector
#' @param method character; "ols" for linear regression or "loclin" for
#'     local linear regression
#' @param w numeric vector; nx1 observation weights,
#'     not used for method "loclin", optional for "ols"
#' @param kernel character; kernel used if method is "loclin". Defaults to
#'     "gauss", the Gaussian kernel.
#' @export
linsmooth = function(X = NULL, y=NULL, method="loclin", w=NULL, kernel="gauss") {
  if(is.null(X) | is.null(y)){
    stop("Must enter values for X or Y")
  }
  if(is.data.frame(X)){
    X = as.matrix(X)
  }
  # If there is no constant column, add one
  colsd = apply(X, 2, var)
  if(all(colsd)!=0){
    X = cbind(rep(1, nrow(X)),X)
  }

  # Perform smoothing methods case-by-case
  if(method=="loclin"){
    if(ncol(X)==2){
      opth = optim(1, function(x) loocv(X, matrix(x), y),
                   lower=0, upper=, method="Brent")$par
    }else{
      opth = optim(1, function(x) loocv(X, diag(x, ncol(X), ncol(X))),
                   method="Nelder-Mead")
    }
    output = loclin_gauss(X, diag(opth, ncol(X), ncol(X)), y)
    return()
  }else if(method=="ols"){
    if(is.null(w)){
      w = rep(1, length(y))
    }
    output = fastols(X, y, w)
  }
  return(output)
}
