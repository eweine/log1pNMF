init_factor_model_log1p <- function(n, p, K) {

  fit <- list()

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

  return(fit)

}
