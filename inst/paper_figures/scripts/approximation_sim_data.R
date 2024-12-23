set.seed(1)
n <- 500
p <- 250
K <- 5

library(distr)
library(Matrix)

l_dist <- UnivarMixingDistribution(
  Unif(0,0.1),
  Exp(rate = 0.5),
  mixCoeff = c(0.9, 0.1)
)

f_dist <- UnivarMixingDistribution(
  Unif(0,0.1),
  Exp(rate = 0.5),
  mixCoeff = c(0.9, 0.1)
)

l_sampler <- distr::r(l_dist)
f_sampler <- distr::r(f_dist)

LL <- matrix(
  data = l_sampler(n = n * K),
  nrow = n,
  ncol = K
)

FF <- matrix(
  data = f_sampler(n = p * K),
  nrow = p,
  ncol = K
)

Lambda <- LL %*% t(FF)

Y <- matrix(
  data = rpois(n = prod(dim(Lambda)), lambda = as.vector(Lambda)),
  nrow = n,
  ncol = p
)

Y_dense <- Y
Y <- as(Y, "CsparseMatrix")

c_vec <- c(
  0.01, 0.1, 0.25, 0.5, 0.75, 1, 2.5,
  5, 7.5, 10, 100, 1000, 10000
)

for (i in 1:length(c_vec)) {

  set.seed(1)
  log1p <- passPCA::fit_factor_model_log1p_exact(
    Y = Y,
    K = 5,
    maxiter = 10000,
    s = rep(c_vec[i], n)
  )

  readr::write_rds(
    log1p,
    glue::glue(
      "~/Documents/data/passPCA/comp_approx/log1p_exact_c{c_vec[i]}.rds"
    )
  )

  set.seed(1)
  log1p_approx <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
    Y = Y,
    K = 5,
    maxiter = 10000,
    s = rep(c_vec[i], n),
    approx_method = "taylor"
  )

  readr::write_rds(
    log1p_approx,
    glue::glue(
      "~/Documents/data/passPCA/comp_approx/log1p_approx_taylor_c{c_vec[i]}.rds"
    )
  )

  Y_tilde <- MatrixExtra::mapSparse(Y, function(x){x/c_vec[i]})
  frob_approx <- NNLM::nnmf(A = as.matrix(Y_tilde), k = 5, max.iter = 10000)

  readr::write_rds(
    frob_approx,
    glue::glue(
      "~/Documents/data/passPCA/comp_approx/frob_approx_c{c_vec[i]}.rds"
    )
  )

}


ll_vec <- numeric(length(c_vec))
ll_vec_log1p_approx <- numeric(length(c_vec))
ll_vec_frobenius_approx <- numeric(length(c_vec))

rmse_vec <- numeric(length(c_vec))
rmse_vec_log1p_approx <- numeric(length(c_vec))
rmse_vec_frobenius_approx <- numeric(length(c_vec))

for (i in 1:length(c_vec)) {

  set.seed(1)

  log1p <- readr::read_rds(
    glue::glue(
      "~/Documents/data/passPCA/comp_approx/log1p_exact_c{c_vec[i]}.rds"
    )
  )

  log1p_approx <- readr::read_rds(
    glue::glue(
      "~/Documents/data/passPCA/comp_approx/log1p_approx_taylor_c{c_vec[i]}.rds"
    )
  )

  frob_approx <- readr::read_rds(
    glue::glue(
      "~/Documents/data/passPCA/comp_approx/frob_approx_c{c_vec[i]}.rds"
    )
  )

  H <- exp(log1p$U %*% t(log1p$V)) - 1
  H_log1p_approx <- exp(log1p_approx$U %*% t(log1p_approx$V)) - 1
  H_frob_approx <- exp(frob_approx$W %*% frob_approx$H) - 1

  ll_vec[i] <- sum(
    dpois(
      x = as.vector(Y_dense),
      lambda = c_vec[i] * as.vector(H),
      log = TRUE
    )
  )

  ll_vec_log1p_approx[i] <- sum(
    dpois(
      x = as.vector(Y_dense),
      lambda = c_vec[i] * as.vector(H_log1p_approx),
      log = TRUE
    )
  )

  ll_vec_frobenius_approx[i] <- sum(
    dpois(
      x = as.vector(Y_dense),
      lambda = c_vec[i] * as.vector(H_frob_approx),
      log = TRUE
    )
  )

  H_df <- data.frame(
    log1p_exact = as.vector(c_vec[i] * as.vector(H)),
    log1p_approx = as.vector(c_vec[i] * as.vector(H_log1p_approx)),
    frob_approx = as.vector(c_vec[i] * as.vector(H_frob_approx))
  )

  rmse_vec[i] <- sqrt(mean((c_vec[i] * as.vector(H) - as.vector(Lambda))^2))
  rmse_vec_log1p_approx[i] <- sqrt(mean((c_vec[i] * as.vector(H_log1p_approx) - as.vector(Lambda))^2))
  rmse_vec_frobenius_approx[i] <- sqrt(mean((c_vec[i] * as.vector(H_frob_approx) - as.vector(Lambda))^2))

}

