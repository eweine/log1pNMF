#' Fit log1p Poisson Factor Model with Exact Log Likelihood
#'
#' @param sc summary of sparse matrix to fit model on
#' @param sc_t summary of transpose of sparse matrix to fit model on
#' @param n number of rows of data matrix
#' @param p number of columns of data matrix
#' @param fit object with model parameters
#' @param maxiter maximum number of updates
#'
#' @return list with LL and FF
#' @keywords internal
#'
#' @useDynLib passPCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
fit_factor_model_log1p_exact <- function(
    sc,
    sc_t,
    n,
    p,
    fit,
    maxiter
) {

  U <- cbind(log(fit$s), fit$LL)
  V <- cbind(rep(1, ncol(Y)), init$FF)
  update_idx <- 1:ncol(fit$LL)

  new_UV <- fit_factor_model_log1p_exact_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    fit$s,
    t(U),
    t(V),
    n,
    p,
    as.integer(maxiter),
    .01,
    .25,
    5,
    update_idx
  )

  fit$LL <- new_UV$U[, -1]
  fit$FF <- new_UV$V[, -1]

  return(fit)

}
