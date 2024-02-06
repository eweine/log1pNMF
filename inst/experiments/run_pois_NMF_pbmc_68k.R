load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/pbmc_68k.RData")
#load("~/Documents/data/fastglmpca/raw_data/pbmc_68k.RData")

set.seed(1)

library(fastTopics)

library(RhpcBLASctl)
omp_set_num_threads(28)
blas_set_num_threads(1)

pois_nmf <- fit_poisson_nmf(
  X = counts,
  k = 10,
  numiter = 1000,
  control = list(nc = 28)
)

readr::write_rds(
  pois_nmf,
  "results/pois_nmf_pbmc_68k_10_factors_1000_iter.rds"
)

