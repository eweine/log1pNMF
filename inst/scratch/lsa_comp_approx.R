# here, I'm interested in investigating the quality of different
# computational approximations to the log1p model.

# I would like to try:
# (1) log1p NMF exact optimization (already done).
# (2) log1p NMF approximate optimization with a Taylor approximation.
# (3) log1p NMF approximate optimization with a Chebyschev approximation.
# (4) vanilla NMF on the shifted log-transformed data.

library(log1pNMF)
library(dplyr)
set.seed(1)
K <- 13

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")
load("../data/experiment_results.Rdata")
barcodes   <- as.data.frame(barcodes)
clusters   <- factor(barcodes$celltype,
                     c("Acinar","Ductal","Endothelial/Mesnchymal","Macrophage",
                       "Alpha","Beta","Delta","Gamma"))

barcodes <- barcodes %>%
  dplyr::mutate(
    condition = if_else(
      condition == "IL-1B_IFNg",
      "IL-1B + IFNg",
      condition
    )
  )

conditions <- factor(barcodes$condition,
                     c("Untreated","IL-1B","IFNg","IL-1B + IFNg"))

log1p_exact <- res_list$pancreas$`1`

############## Chebyschev Approximation ############
set.seed(1)
log1p_cheby <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 13, 
  loglik = "approx",
  approx_technique = "chebyshev",
  control = list(
    maxiter = 100,
    threads = 7
  )
)

set.seed(1)
log1p_cheby_expanded <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 13, 
  loglik = "approx",
  approx_technique = "chebyshev",
  control = list(
    maxiter = 100,
    threads = 7
  ),
  chebyshev_interval = c(.Machine$double.eps, log1p(2.5))
)

set.seed(1)
log1p_cheby_expanded2 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 13, 
  loglik = "approx",
  approx_technique = "chebyshev",
  control = list(
    maxiter = 100,
    threads = 7
  ),
  chebyshev_interval = c(.Machine$double.eps, log1p(1.5))
)

############## NNLM ###########################

s <- Matrix::rowSums(counts)
s <- s / mean(s)
counts <- Matrix::Diagonal(x = (1/s)) %*% counts
counts <- MatrixExtra::mapSparse(counts, log1p)
counts_dense <- as.matrix(counts)

set.seed(1)
nnlm_r1 <- NNLM::nnmf(
  A = counts_dense,
  k = 1
)

L_nnlm_init <- cbind(
  nnlm_r1$W,
  matrix(
    data = 1e-5,
    nrow = nrow(nnlm_r1$W),
    ncol = 12
  )
)

F_nnlm_init <- rbind(
  nnlm_r1$H,
  matrix(
    data = 1e-5,
    ncol = ncol(nnlm_r1$H),
    nrow = 12
  )
)

set.seed(1)
nnlm_full <- NNLM::nnmf(
  A = counts_dense,
  k = 13,
  init = list(W = L_nnlm_init, H = F_nnlm_init)
)

i <- c(sample(which(clusters == "Beta"),900),
       which(clusters != "Beta"))

library(fastTopics)
structure_plot(
  log1pNMF:::normalize_bars(log1p_exact$LL[i,]), 
  grouping = clusters[i],gap = 25,perplexity = 70,n = Inf, font.size = 14
)

structure_plot(
  log1pNMF:::normalize_bars(log1p_cheby$LL[i,]), 
  grouping = clusters[i],gap = 25,perplexity = 70,n = Inf, font.size = 14
)

structure_plot(
  log1pNMF:::normalize_bars(log1p_cheby_expanded$LL[i,]), 
  grouping = clusters[i],gap = 25,perplexity = 70,n = Inf, font.size = 14
)

structure_plot(
  log1pNMF:::normalize_bars(log1p_taylor$LL[i,]), 
  grouping = clusters[i],gap = 25,perplexity = 70,n = Inf, font.size = 14
)

structure_plot(
  log1pNMF:::normalize_bars(nnlm_full$W[i,]), 
  grouping = clusters[i],gap = 25,perplexity = 70,n = Inf, font.size = 14
)

