#' Fit sparse quadratic approximation log1p factor model with C++
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param approx_range range of Chebyschev approximation
#' @param maxiter maximum number of updates
#' @param init_method method for initialization
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
    maxiter = 100,
    init_U = NULL,
    init_V = NULL,
    init_method = c("random", "frob_nmf"),
    s = NULL
) {

  n <- nrow(Y)
  p <- ncol(Y)

  null_s <- FALSE

  if (is.null(s)) {

    s <- rep(1, n)
    null_s <- TRUE

  }

  init_method = match.arg(init_method)

  if (!is.null(init_U) & !is.null(init_V)) {

    init <- list()
    init$U <- init_U
    init$V <- init_V

  } else {

    init <- init_factor_model_log1p(Y, s, n, p, K, init_method)

  }

  sc <- Matrix::summary(Y)
  sc_t <- Matrix::summary(Matrix::t(Y))

  if (!null_s) {

    init$U <- cbind(
      log(s), init$U
    )

    init$V <- cbind(
      rep(1, p), init$V
    )

    update_idx <- 1:K

  } else {

    update_idx <- 0:(K - 1)

  }

  fit <- fit_factor_model_log1p_exact_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    s,
    t(init$U),
    t(init$V),
    n,
    p,
    as.integer(maxiter),
    .01,
    .25,
    5,
    update_idx
  )

  if (!null_s) {

    fit$U <- fit$U[, -1]
    fit$V <- fit$V[, -1]

  }

  return(fit)

}


fit_factor_model_log1p_exact_add_const <- function(
    Y,
    K,
    maxiter,
    cc
) {

  n <- nrow(Y)
  p <- ncol(Y)

  init <- init_factor_model_log1p(Y, s, n, p, K, "random")

  sc <- Matrix::summary(Y)
  sc_t <- Matrix::summary(Matrix::t(Y))

  update_idx <- 0:(K - 1)

  fit <- fit_factor_model_log1p_exact_add_const_cpp_src(
    sc$x,
    sc$i - 1,
    sc$j - 1,
    sc_t$x,
    sc_t$i - 1,
    sc_t$j - 1,
    cc,
    t(init$U),
    t(init$V),
    n,
    p,
    as.integer(maxiter),
    .01,
    .25,
    5,
    update_idx
  )

  return(fit)

}

