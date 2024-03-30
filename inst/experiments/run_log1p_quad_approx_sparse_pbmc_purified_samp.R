load("~/Documents/data/fastglmpca/raw_data/pbmc_purified.RData")
set.seed(1)
counts <- counts[, Matrix::colSums(counts) > 0]

library(passPCA)

set.seed(1)
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 30,
  maxiter = 2,
  approx_range = c(0, 1.25),
  s = as.vector(Matrix::rowSums(counts) / mean(Matrix::rowSums(counts))),
  init_method = "random"
)

library(Matrix)
set.seed(1)
log1p_fit_frob_init <- fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 30,
  maxiter = 2,
  approx_range = c(0, 1.25),
  s = as.vector(Matrix::rowSums(counts) / mean(Matrix::rowSums(counts))),
  init_method = "frob_nmf"
)

mod2 <- readr::read_rds(
  glue::glue(
    "results/log1p_quad_approx_pbmc_purified_30_factors_125_iter.rds"
  )
)

readr::write_rds(
  log1p_fit,
  glue::glue(
    "results/log1p_quad_approx_pbmc_purified_30_factors_125_iter.rds"
    )
)

library(fastTopics)
set.seed(1)
nmf_fit <- fit_poisson_nmf(
  X = counts,
  k = 30,
  numiter = 125
)

