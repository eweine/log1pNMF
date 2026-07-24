#' Fit log1p Poisson Factor Model with Approximate Log Likelihood
#'
#' @param sc_x,sc_i,sc_j nonzero values and 0-based row/column indices of the data
#' @param sc_t_x,sc_t_i,sc_t_j the same, for the transpose of the data
#' @param s size factor
#' @param n number of rows of data matrix
#' @param p number of columns of data matrix
#' @param fit object with model parameters
#' @param maxiter maximum number of updates
#'
#' @return list with fit and progress info
#' @keywords internal
#'
#' @useDynLib log1pNMF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
fit_factor_model_log1p_quad_approx_sparse <- function(
    sc_x,
    sc_i,
    sc_j,
    sc_t_x,
    sc_t_i,
    sc_t_j,
    s,
    n,
    p,
    fit,
    maxiter
) {

  update_idx <- 0:(ncol(fit$LL) - 1)

  # optimization-time budget (seconds); Inf means no limit
  max_time <- if (is.null(fit$control$max_time)) Inf else fit$control$max_time

  new_UV <- fit_factor_model_log1p_quad_approx_sparse_cpp_src(
    sc_x,
    sc_i,
    sc_j,
    sc_t_x,
    sc_t_i,
    sc_t_j,
    s,
    t(fit$LL),
    t(fit$FF),
    fit$a1,
    fit$a2,
    n,
    p,
    as.integer(maxiter),
    fit$control$ls_alpha,
    fit$control$ls_beta,
    fit$control$num_ccd_iter,
    update_idx,
    fit$control$verbose,
    fit$control$tol,
    isTRUE(fit$control$track_time),
    isTRUE(fit$control$track_exact_loglik),
    max_time
  )

  fit$LL <- new_UV$U
  fit$FF <- new_UV$V
  fit$converged <- new_UV$converged
  fit$objective_trace <- new_UV$objective_trace

  if (isTRUE(fit$control$track_time)) {
    fit$time_trace <- new_UV$time_trace
  }

  if (isTRUE(fit$control$track_exact_loglik)) {
    fit$exact_loglik_trace <- new_UV$exact_loglik_trace
  }

  return(fit)

}
