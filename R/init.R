init_factor_model_log1p <- function(Y, s, n, p, K, init_method) {

  fit <- list()

  if (init_method == "random") {

    fit$U <- matrix(
      data = rexp(n * K, rate = 5),
      nrow = n,
      ncol = K
    )

    fit$V <- matrix(
      data = rexp(p * K, rate = 5),
      nrow = p,
      ncol = K
    )

  } else if (init_method == "frob_nmf") {

    Y <- Y / s
    Y <- MatrixExtra::mapSparse(Y, log1p)
    frob_nmf <- RcppML::nmf(
      A = Y, k = K
    )

    d <- sqrt(frob_nmf$d)
    fit$U <- pmin(frob_nmf$w * d + 1e-12, 4)
    fit$V <- pmin(t(d * frob_nmf$h) + 1e-12, 4)

  }

  return(fit)

}
