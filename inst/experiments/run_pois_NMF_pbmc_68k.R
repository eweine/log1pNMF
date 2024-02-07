#load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/pbmc_68k.RData")
load("~/Documents/data/fastglmpca/raw_data/pbmc_68k.RData")

set.seed(1)

library(fastTopics)

pois_nmf <- fit_poisson_nmf(
  X = counts,
  k = 10,
  numiter = 1000,
  control = list(nc = 8)
)

readr::write_rds(
  pois_nmf,
  "results/pois_nmf_pbmc_68k_10_factors_1000_iter.rds"
)

