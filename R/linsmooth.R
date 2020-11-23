#' This function performs linear smoothing
#'
#' @param X numeric matrix; an nxk data matrix of observations
#' @param y numeric vector;an nx1 response vector
#' @param Xeval numeric matrix; mxk data matrix where the linear smoother should
#'     produced predicted values of y. Defaults to X.
#' @param method character; "ols" for linear regression or "loclin" for
#'     local linear regression
#' @param w numeric vector; nx1 observation weights,
#'     not used for method "loclin", optional for "ols"
#' @param kernel character; kernel used if method is "loclin". Defaults to
#'     "gauss", the Gaussian kernel.
#' @param H matrix; Bandwidth matrix used if "loclin" method is used. Defaults
#'     to a diagonal matrix being optimized by minimizing LOOCV score.
#' @export
linsmooth = function(X,
                     y,
                     Xeval = X,
                     method="loclin",
                     w=NULL,
                     kernel="gauss",
                     H=NULL) {
  if(is.data.frame(X)){
    X = as.matrix(X)
  }
  # If there is no constant column, add one
  colsd = apply(X, 2, var)
  if(all(colsd!=0)){
    X = cbind(rep(1, nrow(X)),X)
  }
  k = ncol(X)
  # Perform smoothing methods case-by-case
  if(method=="loclin"){
    if(kernel=="gauss"){
      if(is.null(H)){
        if(ncol(X)==2){
          opth = optim(1, function(h) loocv_gauss(X, matrix(h), y),
                       lower=0, upper=max(X[,2])-min(X[,2]), method="Brent")$par
        }else{
          opth = optim(1, function(h) loocv_gauss(X,diag(h, k-1, k-1), y),
                       method="Nelder-Mead")
        }
      }else{
        opth = H
      }
      output = loclin_gauss(X, matrix(opth), y, Xeval, 1)
      output["bandwidth"] = opth
    }
  }else if(method=="ols"){
    if(is.null(w)){
      w = rep(1, length(y))
    }
    output = fastols(X, y, w)
  }
  return(output)
}