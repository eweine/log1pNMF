#' Fit log1p Poisson Factor Model with Exact Log Likelihood
#'
#' @param sc_x,sc_i,sc_j nonzero values and 0-based row/column indices of the data
#' @param sc_t_x,sc_t_i,sc_t_j the same, for the transpose of the data
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

  # add size factor to L and F
  U <- cbind(log(s), fit$LL)
  V <- cbind(rep(1, p), fit$FF)
  update_idx <- 1:ncol(fit$LL)

  # optimization-time budget (seconds); Inf means no limit
  max_time <- if (is.null(fit$control$max_time)) Inf else fit$control$max_time

  new_UV <- fit_factor_model_log1p_exact_cpp_src(
    sc_x,
    sc_i,
    sc_j,
    sc_t_x,
    sc_t_i,
    sc_t_j,
    s,
    t(U),
    t(V),
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
    max_time
  )

  # remove size factor
  fit$LL <- new_UV$U[, -1, drop = FALSE]
  fit$FF <- new_UV$V[, -1, drop = FALSE]
  fit$converged <- new_UV$converged
  fit$objective_trace <- new_UV$objective_trace

  if (isTRUE(fit$control$track_time)) {
    fit$time_trace <- new_UV$time_trace
  }

  # the exact objective already is the model log-likelihood core; expose it
  # under exact_loglik_trace when requested so downstream code is uniform
  # across the approximate and exact fitting paths
  if (isTRUE(fit$control$track_exact_loglik)) {
    fit$exact_loglik_trace <- new_UV$objective_trace
  }

  return(fit)

}
