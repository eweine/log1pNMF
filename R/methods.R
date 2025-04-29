#' @rdname fitted.log1p_nmf_fit
#'
#' @title Get Fitted Values (i.e. values of \eqn{\lambda}) for log1p NMF Model
#'
#' @description \code{fitted} method for the
#'   \dQuote{log1p_nmf} class.
#'
#' @param object An object of class \dQuote{log1p_nmf_fit},
#'   typically the result of calling \code{\link{fit_poisson_log1p_nmf}}.
#'
#' @param \dots Additional arguments passed to the generic
#'   \code{fitted} method (that will not do anything).
#'
#' @return An n x p matrix of fitted means. Calculated as
#'   \deqn{\exp(LF')} using the \code{fit} object.
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
  return(Lambda)

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
#' @return A double indicating the log-likelihood of the fitted model.
#'
#' @method logLik log1p_nmf_fit
#'
#' @export
#'
logLik.log1p_nmf_fit <- function (object, Y, ...) {
  if (!inherits(object, "log1p_nmf_fit")) {
    
    stop("object must be of class log1p_nmf_fit")
    
  }
  
  if (!inherits(Y, "dgCMatrix")) {
    
    Y <- as(Y, "CsparseMatrix")
    
  }
  
  n <- nrow(Y)
  p <- ncol(Y)
  sqrt_alpha <- sqrt(max(1, object$cc))
  s <- object$s * object$cc
  
  U <- cbind(log(s), object$LL)
  V <- cbind(rep(1, p), object$FF)

  sc <- Matrix::summary(Y)
  ll <- get_loglik_exact(
    sqrt_alpha * t(U),
    sqrt_alpha * t(V),
    sc$x,
    sc$i - 1,
    sc$j - 1,
    s,
    n,
    p
  ) + sum(s) - sum(lfactorial(Y@x))
  
  return(ll)
  
}
