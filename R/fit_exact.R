#' Fit sparse quadratic approximation log1p factor model with C++
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param cc value of c.
#' @param fit_c if c should be fit via maximum likelihood.
#' @param approx_range range of Chebyschev approximation
#' @param maxiter maximum number of updates
#'
#' @return list with fit and progress info
#' @export
#'
#' @useDynLib passPCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
fit_factor_model_log1p_exact <- function(
    Y,
    K,
    maxiter,
    s = NULL,
    cc = 1,
    fit_c = FALSE
) {

  n <- nrow(Y)
  p <- ncol(Y)

  if (is.null(s)) {

    s <- rep(1, n)

  }

  init <- init_factor_model_log1p(n, p, K)

  init$U <- cbind(
    log(s), init$U
  )

  init$V <- cbind(
    rep(1, p), init$V
  )

  sc <- Matrix::summary(Y)
  sc_t <- Matrix::summary(Matrix::t(Y))

  # pre-compute constants for fitting c
  sum_Y <- sum(sc$x)
  sum_s <- sum(s)

  fit <- fit_factor_model_log1p_exact_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    sum_Y,
    sum_s,
    s,
    cc,
    t(init$U),
    t(init$V),
    n,
    p,
    as.integer(maxiter),
    .01,
    .25,
    5,
    1:K,
    fit_c
  )

  # remove size factor after fitting
  fit$U <- fit$U[, -1]
  fit$V <- fit$V[, -1]

  return(fit)

}
