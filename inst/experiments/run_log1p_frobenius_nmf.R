#load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")
load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

counts <- Matrix::rowScale(counts, 1 / Matrix::rowSums(counts))
counts <- MatrixExtra::mapSparse(counts, log1p)

library(RcppML)

nmf_out <- nmf(counts, k = 10)
readr::write_rds(
  nmf_out,
  "~/Documents/passPCA/inst/experiments/results/log1p_frobenius_nmf_droplets.rds"
)

