#' This function performs linear smoothing
#'
#' @param X numeric matrix; an nxk data matrix of observations
#' @param y numeric vector;an nx1 response vector
#' @param Xeval numeric matrix; mxk data matrix where the linear smoother should
#'     produced predicted values of y. Defaults to X.
#' @param method character; "ols" for linear regression or "loclin" for
#'     local linear regression. Defaults to "ols"
#' @param w numeric vector; nx1 observation weights,
#'     not used for method "loclin", optional for "ols"
#' @param kernel character; kernel used if method is "loclin". Defaults to
#'     "gauss", the Gaussian kernel.
#' @param H positive vector; Bandwidth vector that represents diagonal of
#'     bandwidth matrix used if "loclin" method is used. Defaults
#'     to a diagonal matrix being optimized by minimizing LOOCV score.
#' @param g an nx1 vector of groups if performing linear smoothing within
#'     groups.
#' @param compute_se logical; If TRUE then OLS procedure will compute standard
#'     errors. TRUE by default.
#' @param compute_hat logical; If TRUE then diagonal of hat matrix (leverage) is
#'     calculated for each observation. FALSE by default.
#' @param nthreads positive integer; Number of threads to use for parallelized
#'     operations. Defaults to 1 thread.
#' @export
linsmooth = function(X,
                     y,
                     Xeval=NULL,
                     method="ols",
                     w=NULL,
                     kernel="gauss",
                     H=NULL,
                     bygroup=NULL,
                     compute_se=TRUE,
                     compute_hat=FALSE,
                     nthreads=1L) {
  # Control input ==============================================================
  if(!is.matrix(X)){
    X = unname(as.matrix(X))
  }
  k = ncol(X)
  n = nrow(X)
  if(!is.numeric(y)){
    stop("Input y must be a numeric vector.")
  }
  if(length(y)!=n){
    stop("Input y must have the same number of elements as the number of rows of X.")
  }
  if(is.null(Xeval)){
    Xevalsame = TRUE
  }else {
    if(!is.matrix(Xeval)){
      Xeval = unname(as.matrix(Xeval))
    }
    if(ncol(Xeval)!=k){
      stop("Xeval must have the same number of columns as X.")
    }
    Xevalsame = FALSE
  }

  compute_se  = as.integer(compute_se)
  compute_hat = as.integer(compute_hat)

  # If there is no constant column, add one
  colsd = apply(X, 2, var)
  if(all(colsd!=0)){
    X = cbind(rep(1, n), X)
    if(!Xevalsame){
      Xeval = cbind(rep(1, nrow(Xeval)), Xeval)
    }
    k = k+1
  }

  # Weights
  if(!is.numeric(w)){
    stop("w must be a vector of positive numeric weights.")
  }

  # Perform smoothing methods case-by-case =====================================
  if(method=="loclin"){
    kern = ifelse(kernel=="gauss", 1L, 2L)
    if(kern!=1L & k>2){
      stop("Multivariate smoothing only allowed for kernal='gauss'")
    }
    if(is.null(H)){
      if(ncol(X)==2){
        H = optim(1, function(h) loclin(X, matrix(h), y, X, sameX=1, kernel=kern)[["cvscore"]],
                     lower=0, upper=max(X[,2])-min(X[,2]), method="Brent", control=list(maxit=8))$par
      }else{
        H = optim(1, function(h) loclin(X,diag(h, k-1, k-1), y, X, sameX=1, kernel=kern)[["cvscore"]],
                     method="Nelder-Mead", control=list(maxit=10))$par
      }
    }
    output = loclin(X, diag(H, k-1, k-1), y, Xeval, 1, kern)
    output["bandwidth"] = H
    return(output)
  }else if(method=="ols"){
    if(is.null(w)){
      w = rep(1, length(y))
    }
    if(!is.null(bygroup)){
      ind = order(g)
      gsorted = as.integer(cumsum(!duplicated(g[ind])))
      output = fastols_by(X[ind,], y[ind], w[ind], gsorted, nthr=nthreads,
                          compute_se = compute_se, compute_hat = compute_hat)
    }else{
      fastols(x,y,w)
    }
    return(output)
  }
}
