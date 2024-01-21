#' Fit sparse linear approximation log1p factor model
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param a slope of linear approximation
#' @param maxiter maximum number of updates
#'
#' @return list with fit and progress info
#' @export
#'
#' @useDynLib passPCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
fit_factor_model_log1p_lin_approx_sparse <- function(
    Y,
    K,
    a = 1,
    maxiter = 25
) {

  n <- nrow(Y)
  p <- ncol(Y)

  fit <- init_factor_model_log1p(n, p, K)

  count_data <- prep_data(Y)

  loglik <- get_loglik_lin_approx_sparse(
    fit$U,
    fit$V,
    count_data$sc$x,
    count_data$sc$i - 1,
    count_data$sc$j - 1,
    a
  )

  iter <- 0

  cat(sprintf("Fitting log1p factor model to %d x %d count matrix.\n",n,p))

  loglik_history <- numeric(maxiter + 1)
  loglik_history[1] <- loglik

  for (iter in 1:maxiter) {

    cat(sprintf("Iteration %d: objective = %+0.12e\n",
                iter - 1, loglik))

    new_U_T <- regress_cols_of_Y_on_X_log1p_lin_approx_sparse(
      t(fit$V),
      count_data$y_rows_data,
      count_data$y_rows_idx,
      a * colSums(fit$V),
      a,
      t(fit$U),
      0:(K - 1),
      5,
      .01,
      .25
    )

    fit$U <- t(new_U_T)

    new_V_T <- regress_cols_of_Y_on_X_log1p_lin_approx_sparse(
      t(fit$U),
      count_data$y_cols_data,
      count_data$y_cols_idx,
      a * colSums(fit$U),
      a,
      t(fit$V),
      0:(K - 1),
      5,
      .01,
      .25
    )

    fit$V <- t(new_V_T)

    loglik <- get_loglik_lin_approx_sparse(
      fit$U,
      fit$V,
      count_data$sc$x,
      count_data$sc$i - 1,
      count_data$sc$j - 1,
      a
    )

    loglik_history[iter + 1] <- loglik

  }

  fit$loglik <- loglik_history
  return(fit)

}
