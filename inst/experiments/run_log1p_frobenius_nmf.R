#load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")
load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

counts <- Matrix::rowScale(counts, 1 / Matrix::rowSums(counts))
counts <- MatrixExtra::mapSparse(counts, log1p)

library(Matrix)
library(RcppML)

nmf_out <- nmf(counts, k = 7)
readr::write_rds(
  nmf_out,
  "~/Documents/passPCA/inst/experiments/results/log1p_frobenius_nmf_droplets.rds"
)

LL <- nmf_out$w
FF <- t(diag(nmf_out$d) %*% nmf_out$h)

rownames(FF) <- colnames(counts)

LL_max <- apply(LL,2,max)

LL_new <- LL %*% diag(1/LL_max)
FF_new <- FF %*% diag(LL_max)

g1 <- structure_plot(LL_new, grouping = samples$tissue)

# now, read in vanilla model
pois_nmf <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/")
