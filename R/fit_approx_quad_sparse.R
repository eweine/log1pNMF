#' Fit sparse quadratic approximation log1p factor model
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param approx_range range of Chebyschev approximation
#' @param maxiter maximum number of updates
#'
#' @return list with fit and progress info
#' @export
#'
fit_factor_model_log1p_quad_approx_sparse <- function(
    Y,
    K,
    approx_range,
    maxiter
) {

  n <- nrow(Y)
  p <- ncol(Y)

  fit <- init_factor_model_log1p(n, p, K)

  count_data <- prep_data(Y)

  # get the approximation
  poly_approx <- pracma::polyApprox(
    exp,
    approx_range[1],
    approx_range[2],
    2
  )

  a1 <- poly_approx$p[2]
  a2 <- poly_approx$p[1]

  loglik <- get_loglik_quad_approx_sparse(
    fit$U,
    fit$V,
    count_data$sc$x,
    count_data$sc$i - 1,
    count_data$sc$j - 1,
    a1,
    a2
  )

  loglik_exact <- get_loglik_exact(
    t(fit$U),
    t(fit$V),
    count_data$sc$x,
    count_data$sc$i - 1,
    count_data$sc$j - 1,
    n,
    p
  )

  iter <- 0

  cat(sprintf("Fitting log1p factor model to %d x %d count matrix.\n",n,p))

  loglik_history <- numeric(maxiter + 1)
  loglik_history[1] <- loglik

  loglik_exact_history <- numeric(maxiter + 1)
  loglik_exact_history[1] <- loglik_exact

  for (iter in 1:maxiter) {

    cat(sprintf("Iteration %d: objective = %+0.12e\n",
                iter - 1, loglik))

    new_U_T <- regress_cols_of_Y_on_X_log1p_quad_approx_sparse(
      t(fit$V),
      count_data$y_rows_data,
      count_data$y_rows_idx,
      a1 * colSums(fit$V),
      crossprod(fit$V),
      t(fit$U),
      a1,
      a2,
      0:(K - 1),
      5,
      .01,
      .25
    )

    fit$U <- t(new_U_T)

    new_V_T <- regress_cols_of_Y_on_X_log1p_quad_approx_sparse(
      t(fit$U),
      count_data$y_cols_data,
      count_data$y_cols_idx,
      a1 * colSums(fit$U),
      crossprod(fit$U),
      t(fit$V),
      a1,
      a2,
      0:(K - 1),
      5,
      .01,
      .25
    )

    fit$V <- t(new_V_T)

    loglik <- get_loglik_quad_approx_sparse(
      fit$U,
      fit$V,
      count_data$sc$x,
      count_data$sc$i - 1,
      count_data$sc$j - 1,
      a1,
      a2
    )

    loglik_exact <- get_loglik_exact(
      t(fit$U),
      t(fit$V),
      count_data$sc$x,
      count_data$sc$i - 1,
      count_data$sc$j - 1,
      n,
      p
    )

    loglik_history[iter + 1] <- loglik
    loglik_exact_history[iter + 1] <- loglik_exact

  }

  fit$loglik <- loglik_history
  fit$loglik_exact <- loglik_exact_history
  return(fit)

}
