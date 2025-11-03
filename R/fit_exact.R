#' Fit log1p Poisson Factor Model with Exact Log Likelihood
#'
#' @param sc summary of sparse matrix to fit model on
#' @param sc_t summary of transpose of sparse matrix to fit model on
#' @param s size factor
#' @param n number of rows of data matrix
#' @param p number of columns of data matrix
#' @param fit object with model parameters
#' @param maxiter maximum number of updates
#'
#' @return list with LL and FF
#' @keywords internal
#'
#' @useDynLib log1pNMF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
fit_factor_model_log1p_exact <- function(
    sc,
    sc_t,
    s,
    n,
    p,
    fit,
    maxiter
) {

  # add size factor to L and F
  U <- cbind(log(s), fit$LL)
  V <- cbind(rep(1, p), fit$FF)

  new_UV <- fit_factor_model_log1p_exact_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    s,
    t(U),
    t(V),
    n,
    p,
    as.integer(maxiter),
    fit$control$ls_alpha,
    fit$control$ls_beta,
    fit$control$num_ccd_iter,
    fit$update_idx_LL,
    fit$update_idx_FF,
    fit$control$verbose,
    fit$control$tol
  )

  # remove size factor
  fit$LL <- new_UV$U[, -1, drop = FALSE]
  fit$FF <- new_UV$V[, -1, drop = FALSE]
  fit$converged <- new_UV$converged
  fit$objective_trace <- new_UV$objective_trace

  return(fit)

}
