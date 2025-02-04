#' Fit log1p Poisson Factor Model with Approximate Log Likelihood
#'
#' @param sc summary of sparse matrix to fit model on
#' @param sc_t summary of transpose of sparse matrix to fit model on
#' @param n number of rows of data matrix
#' @param p number of columns of data matrix
#' @param fit object with model parameters
#' @param maxiter maximum number of updates
#'
#' @return list with fit and progress info
#' @keywords internal
#'
#' @useDynLib passPCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
fit_factor_model_log1p_quad_approx_sparse <- function(
    sc,
    sc_t,
    n,
    p,
    fit,
    maxiter
) {

  update_idx <- 0:(ncol(fit$LL) - 1)

  new_UV <- fit_factor_model_log1p_quad_approx_sparse_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    s,
    t(fit$LL),
    t(fit$FF),
    fit$a1,
    fit$a2,
    n,
    p,
    as.integer(maxiter),
    .01,
    .25,
    5,
    update_idx
  )

  fit$LL <- new_UV$U
  fit$FF <- new_UV$V

  return(fit)

}
