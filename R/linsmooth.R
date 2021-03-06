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
#' @param by an nx1 vector of groups if performing linear smoothing within
#'     groups.
#' @param compute_se logical; If TRUE then OLS procedure will compute standard
#'     errors. TRUE by default.
#' @param compute_hat logical; If TRUE then diagonal of hat matrix (leverage) is
#'     calculated for each observation. FALSE by default.
#' @param nthr positive integer; Number of threads to use for parallelized
#'     operations. Defaults to 1 thread.
#' @export
linsmooth = function(X,
                     y,
                     Xeval=NULL,
                     method="ols",
                     w=NULL,
                     kernel="gauss",
                     H=NULL,
                     by=NULL,
                     compute_se=TRUE,
                     compute_hat=FALSE,
                     nthr=1L) {
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
  colsd = apply(X, 2, stats::var)
  if(all(colsd!=0)){
    X = cbind(rep(1, n), X)
    if(!Xevalsame){
      Xeval = cbind(rep(1, nrow(Xeval)), Xeval)
    }
    k = k+1
  }

  # Weights
  if(!is.null(w) & !is.numeric(w)){
    stop("w must be a vector of positive numeric weights.")
  }

  # Perform smoothing methods case-by-case =====================================
  if(method=="loclin"){
    kern = ifelse(kernel=="gauss", 1L,
                  ifelse(kernel=="epan", 2L, 3L))
    if(kern!=1L & k>2){
      stop("Multivariate smoothing only allowed for kernal='gauss'")
    }
    if(is.null(H)){
      if(ncol(X)==2){
        H = stats::optim(1, function(h) loclin_sameX(X, y, matrix(h), kernel=kern)[["cvscore"]],
                         lower=0, upper=max(X[,2])-min(X[,2]), method="Brent", control=list(maxit=8))$par
      }else{
        H = stats::optim(1, function(h) loclin_sameX(X, y, diag(h, k-1, k-1), kernel=kern)[["cvscore"]],
                         method="Nelder-Mead", control=list(maxit=10))$par
      }
    }
    if(is.null(Xeval)){
      output = loclin_sameX(X, y, diag(H, k-1, k-1), kern)
    }else {
      output = loclin_sameX(X, y, diag(H, k-1, k-1), Xeval, kern)
    }
    output["bandwidth"] = H
    return(output)
  }else if(method=="ols"){
    if(!is.null(by)){
      ind = order(by)
      gsorted = as.integer(cumsum(!duplicated(by[ind])))
      if(!is.null(w)){
        output = fastols_bywt(X[ind,], y[ind], w[ind], gsorted, nthr=nthr,
                              compute_se = compute_se, compute_hat = compute_hat)
      }else{
        output = fastols_by(X[ind,], y[ind], gsorted, nthr=nthr,
                            compute_se = compute_se, compute_hat = compute_hat)
      }
    }else{
      if(!is.null(w)){
        output = fastolswt(X, y, w,
                           compute_se = compute_se, compute_hat = compute_hat)
      }else{
        output = fastols(X, y,
                         compute_se = compute_se, compute_hat = compute_hat)
      }
    }
    return(output)
  }
}

