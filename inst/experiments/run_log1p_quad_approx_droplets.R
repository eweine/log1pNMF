load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")
#load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

library(RhpcBLASctl)
omp_set_num_threads(28)
blas_set_num_threads(1)

set.seed(1)

library(passPCA)

log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 10,
  maxiter = 5,
  approx_range = c(0, 1.25),
  s = Matrix::rowSums(counts) / mean(Matrix::rowSums(counts))
)

readr::write_rds(
  log1p_fit,
  "results/log1p_quad_approx_droplets_10_factors_2500_iter.rds"
)
