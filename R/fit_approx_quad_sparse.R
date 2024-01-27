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
fit_factor_model_log1p_quad_approx_sparse <- function(
    Y,
    K,
    approx_range,
    maxiter
) {

  n <- nrow(Y)
  p <- ncol(Y)

  init <- init_factor_model_log1p(n, p, K)

  # get the approximation
  poly_approx <- pracma::polyApprox(
    exp,
    approx_range[1],
    approx_range[2],
    2
  )

  a1 <- poly_approx$p[2]
  a2 <- poly_approx$p[1]

  sc <- Matrix::summary(Y)
  sc_t <- Matrix::summary(Matrix::t(Y))

  fit <- fit_factor_model_log1p_quad_approx_sparse_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    t(init$U),
    t(init$V),
    a1,
    a2,
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
