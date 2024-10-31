#' Fit sparse quadratic approximation log1p factor model with C++
#'
#' @param Y sparse matrix.
#' @param K rank of factorization
#' @param approx_range range of Chebyschev approximation
#' @param maxiter maximum number of updates
#' @param init_method method for initialization
#' @param fit_constant boolean indicating if constant should be fit
#' @param constant_init initialization value for constant
#' @param tol tolerance for optimization
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
    init_method = c("random", "frob_nmf"),
    s = NULL,
    fit_constant = FALSE,
    constant_init = NULL,
    tol = 1e-4
) {

  n <- nrow(Y)
  p <- ncol(Y)

  null_s <- FALSE

  if (is.null(s)) {

    s <- rep(1, n)
    null_s <- TRUE

  }

  init_method = match.arg(init_method)
  init <- init_factor_model_log1p(Y, s, n, p, K, init_method)

  sc <- Matrix::summary(Y)
  sc_t <- Matrix::summary(Matrix::t(Y))

  if (!null_s) {

    init$U <- cbind(
      log(s), init$U
    )

    init$V <- cbind(
      rep(1, p), init$V
    )

  }

  if (fit_constant) {

    if (!is.null(constant_init)) {

      init$U <- cbind(
        rep(constant_init, n), init$U
      )

    } else {

      init$U <- cbind(
        rep(0, n), init$U
      )

    }

    init$V <- cbind(
      rep(1, p), init$V
    )

  } else if (!is.null(constant_init)) {

    init$U <- cbind(
      rep(constant_init, n), init$U
    )

    init$V <- cbind(
      rep(1, p), init$V
    )

  }

  K_total <- ncol(init$U)
  K_not_fit <- K_total - K

  update_idx <- K_not_fit:(K_total - 1)

  print("update_idx = ")
  print(update_idx)

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
    update_idx,
    fit_constant,
    tol
  )

  # here, I need to adjust what I'm doing in order to extract the right fit
  # I have to be a bit careful about this but I think it should be fairly easy

  #fit$U <- fit$U[,(update_idx + 1),drop = FALSE]
  #fit$V <- fit$V[,(update_idx + 1),drop = FALSE]

  return(fit)

}
