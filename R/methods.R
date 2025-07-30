#' @rdname fitted.log1p_nmf_fit
#'
#' @title Get Fitted Values for log1p NMF Model
#'
#' @description \code{fitted} method for the
#'   \code{log1p_nmf} class. Computes a matrix of values where entry i, j is
#' \eqn{s_{i}\lambda_{ij}}.
#'
#' @param object An object of class \dQuote{log1p_nmf_fit},
#'   typically the result of calling \code{\link{fit_poisson_log1p_nmf}}.
#'
#' @param \dots Additional arguments passed to the generic
#'   \code{fitted} method (that will not do anything).
#'
#' @return An n x p matrix of fitted means. 
#'
#' @method fitted log1p_nmf_fit
#'
#' @export
#'
fitted.log1p_nmf_fit <- function (object, ...) {
  if (!inherits(object, "log1p_nmf_fit")) {

    stop("object must be of class log1p_nmf_fit")

  }

  Lambda <- object$cc * (exp(tcrossprod(object$LL, object$FF) / object$alpha) - 1)
  fitted_vals <- object$s * Lambda
  return(fitted_vals)

}

#' @rdname logLik.log1p_nmf_fit
#'
#' @title Get log likelihood of data Y for fitted log1p NMF Model
#'
#' @description \code{logLik} method for the
#'   \dQuote{log1p_nmf} class.
#'
#' @param object An object of class \dQuote{log1p_nmf_fit},
#'   typically the result of calling \code{\link{fit_poisson_log1p_nmf}}.
#'
#' @param Y data matrix. Must be convertible to an object of
#' type \code{dgCMatrix}.
#' 
#' @param s custom size factor. Should be provided if calculating the
#' log-likelihood on data that the model was not trained on.
#'
#' @return A double indicating the log-likelihood of the fitted model.
#'
#' @method logLik log1p_nmf_fit
#'
#' @export
#'
logLik.log1p_nmf_fit <- function (object, Y, s = NULL, ...) {
  if (!inherits(object, "log1p_nmf_fit")) {
    
    stop("object must be of class log1p_nmf_fit")
    
  }
  
  if (!inherits(Y, "dgCMatrix")) {
    
    Y <- as(Y, "CsparseMatrix")
    
  }
  
  n <- nrow(Y)
  p <- ncol(Y)
  sqrt_alpha <- sqrt(max(1, object$cc))
  
  if (!is.null(s)) {
    
    if (!is.positive.numeric(s) || length(s) != nrow(Y)) {
      
      stop("s must be a positive vector of length nrow(Y)")
      
    }
    
  } else {
    
    s <- object$s
    
  }
  
  s <- s * object$cc
  
  U <- cbind(log(s), (1 / sqrt_alpha) * object$LL)
  V <- cbind(rep(1, p), (1 / sqrt_alpha) * object$FF)

  sc <- Matrix::summary(Y)
  ll <- get_loglik_exact(
    t(U),
    t(V),
    sc$x,
    sc$i - 1,
    sc$j - 1,
    s,
    n,
    p
  ) + n*p*object$cc - sum(lfactorial(Y@x))
  
  return(ll)
  
}
