init_factor_model_log1p <- function(n, p, K, s) {

  fit <- list()

  fit$U <- cbind(
    log(s),
    matrix(
      data = rexp(n * K, rate = 15),
      nrow = n,
      ncol = K
    )
  )

  fit$V <- cbind(
    rep(1, p),
    matrix(
      data = rexp(p * K, rate = 15),
      nrow = p,
      ncol = K
    )
  )

  return(fit)

}
