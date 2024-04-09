#' Fit sparse quadratic approximation log1p factor model with C++
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param approx_range range of Chebyschev approximation
#' @param maxiter maximum number of updates
#' @param init_U initialization of U
#' @param init_V initialization of V
#' @param update_idx indication of which indices should be updated
#' @param s size factor
#' @param init_method method for initialization
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
    maxiter,
    init_U = NULL,
    init_V = NULL,
    update_idx = NULL,
    init_method = c("random", "frob_nmf"),
    s = NULL
) {

  init_method = match.arg(init_method)

  n <- nrow(Y)
  p <- ncol(Y)

  if (is.null(s)) {

    s <- rep(1, n)

  }

  if (!is.null(init_U) & !is.null(init_V)) {

    init <- list()
    init$U <- init_U
    init$V <- init_V

  } else {

    init <- init_factor_model_log1p(Y, s, n, p, K, init_method)

  }

  if (is.null(update_idx)) {

    C_update_idx <- 0:(K - 1)

  } else {

    C_update_idx <- update_idx - 1

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

  fit <- fit_factor_model_log1p_quad_approx_sparse_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    s,
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
    C_update_idx
  )

  return(fit)

}
