load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")
#load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

set.seed(1)

library(fastTopics)

pois_nmf <- fit_poisson_nmf(
  X = counts,
  k = 10,
  numiter = 2500,
  control = list(nc = 28)
)

readr::write_rds(
  pois_nmf,
  "results/pois_nmf_droplets_10_factors_2500_iter.rds"
)
