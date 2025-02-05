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
#'   \code{fitted} method.
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

  Lambda <- object$cc * (exp(tcrossprod(object$LL, object$FF)) - 1)
  return(Lambda)

}
