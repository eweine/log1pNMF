#' Fit sparse quadratic approximation log1p factor greedily
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param approx_range range of Chebyschev approximation
#' @param iter_per_factor number of updates per factors
#' @param s size factor
#'
#' @return list with fit and progress info
#' @export
#'
#' @useDynLib passPCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
fit_factor_model_log1p_quad_approx_sparse_greedy <- function(
    Y,
    K,
    approx_range,
    iter_per_factor,
    s = NULL
) {

  n <- nrow(Y)
  p <- ncol(Y)

  if (is.null(s)) {

    s <- rep(1, n)

  }

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

  fit <- list()

  fit$U <- matrix(
    nrow = n,
    ncol = 1,
    data = rexp(n, rate = 15)
  )

  fit$V <- matrix(
    nrow = p,
    ncol = 1,
    data = rexp(p, rate = 15)
  )

  for (k in 1:K) {

    fit <- fit_factor_model_log1p_quad_approx_sparse_cpp_src(
      sc$x,
      sc$i - 1,
      sc$j - 1,
      sc_t$x,
      sc_t$i - 1,
      sc_t$j - 1,
      s,
      t(fit$U),
      t(fit$V),
      a1,
      a2,
      n,
      p,
      as.integer(iter_per_factor),
      .01,
      .25,
      5,
      k - 1
    )

    fit$U <- cbind(fit$U, rexp(n, 15))
    fit$V <- cbind(fit$V, rexp(p, 15))

  }

  return(fit)

}
