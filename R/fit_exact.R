#' Fit exact log1p factor model
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param maxiter maximum number of updates
#'
#' @return list with fit and progress info
#' @export
#'
#' @useDynLib passPCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
fit_factor_model_log1p <- function(
  Y,
  K,
  maxiter
) {

  n <- nrow(Y)
  p <- ncol(Y)

  fit <- init_factor_model_log1p(n, p, K)

  count_data <- prep_data(Y)

  loglik <- get_loglik_exact(
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

  for (iter in 1:maxiter) {

    cat(sprintf("Iteration %d: objective = %+0.12e\n",
                iter - 1, loglik))

    new_U_T <- regress_cols_of_Y_on_X_log1p_pois_exact(
      fit$V,
      count_data$y_rows_data,
      count_data$y_rows_idx,
      t(fit$U),
      0:(K - 1),
      5,
      .01,
      .25
    )

    fit$U <- t(new_U_T)

    new_V_T <- regress_cols_of_Y_on_X_log1p_pois_exact(
      fit$U,
      count_data$y_cols_data,
      count_data$y_cols_idx,
      t(fit$V),
      0:(K - 1),
      5,
      .01,
      .25
    )

    fit$V <- t(new_V_T)

    loglik <- get_loglik_exact(
      new_U_T,
      new_V_T,
      count_data$sc$x,
      count_data$sc$i - 1,
      count_data$sc$j - 1,
      n,
      p
    )

    loglik_history[iter + 1] <- loglik

  }

  fit$loglik <- loglik_history
  return(fit)

}


#' Fit sparse quadratic approximation log1p factor model with C++
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param approx_range range of Chebyschev approximation
#' @param maxiter maximum number of updates
#'
#' @return list with fit and progress info
#' @export
#'
#' @useDynLib passPCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
fit_factor_model_log1p_exact_cpp <- function(
    Y,
    K,
    maxiter
) {

  n <- nrow(Y)
  p <- ncol(Y)

  init <- init_factor_model_log1p(n, p, K)

  sc <- Matrix::summary(Y)
  sc_t <- Matrix::summary(Matrix::t(Y))

  fit <- fit_factor_model_log1p_exact_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    t(init$U),
    t(init$V),
    n,
    p,
    as.integer(maxiter),
    .01,
    .25,
    5,
    0:(K-1)
  )

  return(fit)

}
