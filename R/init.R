init_factor_model_log1p <- function(n, p, K) {

  fit <- list()

  fit$U <- matrix(
    data = rexp(n * K, rate = 15),
    nrow = n,
    ncol = K
  )

  fit$V <- matrix(
    data = rexp(p * K, rate = 15),
    nrow = p,
    ncol = K
  )

  return(fit)

}
