fit_nn_glmpca <- function(
    Y,
    K,
    maxiter
) {

  n <- nrow(Y)
  p <- ncol(Y)

  init <- init_factor_model_log1p(Y, s, n, p, K, "random")

  sc <- Matrix::summary(Y)

  update_idx <- 0:(K - 1)

  fit <- fit_factor_model_nn_glmpca_cpp_src(
    Y,
    sc$x,
    sc$i - 1,
    sc$j - 1,
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
