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
fit_factor_model_log1p_lin_approx_sparse <- function(
    Y,
    K,
    a,
    maxiter,
    s = NULL
) {

  n <- nrow(Y)
  p <- ncol(Y)

  if (is.null(s)) {

    s <- rep(1, n)

  }

  init <- init_factor_model_log1p(n, p, K, s)

  sc <- Matrix::summary(Y)
  sc_t <- Matrix::summary(Matrix::t(Y))

  fit <- fit_factor_model_log1p_lin_approx_sparse_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    s,
    t(init$U),
    t(init$V),
    a,
    n,
    p,
    as.integer(maxiter),
    .01,
    .25,
    5,
    1:K
  )

  fit$U <- fit$U[, -1]
  fit$V <- fit$V[, -1]
  fit$s <- s
  fit$loglik <- fit$loglik - sum(MatrixExtra::mapSparse(Y, lfactorial))

  return(fit)

}
