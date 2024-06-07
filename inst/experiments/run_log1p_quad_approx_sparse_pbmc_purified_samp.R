load("~/Documents/data/fastglmpca/raw_data/pbmc_purified.RData")
set.seed(1)
sz <- Matrix::rowSums(counts) / mean(Matrix::rowSums(counts))
counts <- counts[, Matrix::colSums(counts) >= 100]

library(passPCA)

set.seed(1)
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse_greedy(
  Y = counts,
  K = 25,
  iter_per_factor = 8,
  approx_range = c(0, 1.25),
  s = as.vector(sz)
)

rownames(log1p_fit$V) <- colnames(counts)
rownames(log1p_fit$U) <- rownames(counts)

readr::write_rds(log1p_fit, "~/Documents/passPCA/inst/experiments/results/pbmc_purified_25K_greedy_log1p_fit.rds")

log1p_fit_init_greedy <- fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 25,
  maxiter = 50,
  approx_range = c(0, 1.25),
  s = as.vector(sz),
  init_U = log1p_fit$U,
  init_V = log1p_fit$V
)

rownames(log1p_fit_init_greedy$V) <- colnames(counts)
rownames(log1p_fit_init_greedy$U) <- rownames(counts)

readr::write_rds(log1p_fit_init_greedy, "~/Documents/passPCA/inst/experiments/results/pbmc_purified_25K_greedy_init_log1p_fit.rds")

library(Matrix)
set.seed(1)
log1p_fit_init_nmf <- fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 25,
  maxiter = 25,
  approx_range = c(0, 1.25),
  s = as.vector(sz),
  init_method = "frob_nmf"
)

rownames(log1p_fit_init_nmf$V) <- colnames(counts)
rownames(log1p_fit_init_nmf$U) <- rownames(counts)

readr::write_rds(log1p_fit_init_nmf, "~/Documents/passPCA/inst/experiments/results/pbmc_purified_25K_frob_nmf_init_log1p_fit.rds")
