library(Matrix)
library(log1pNMF)
library(dplyr)

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")

library(NNLM)

set.seed(1)

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

i <- c(sample(which(clusters == "Beta"),900),
       which(clusters != "Beta"))

K <- 12

s <- Matrix::rowSums(counts) 
s <- s / mean(s)
Y_tilde <- Matrix::Diagonal(x = 1/s) %*% counts
Y_tilde <- MatrixExtra::mapSparse(Y_tilde, log1p)

set.seed(1)
fit0 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 13,
  init_method = "rank1",
  loglik = "exact",
  control = list(maxiter = 1)
)

fit_approx <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 13,
  loglik = "approx",
  control = list(maxiter = 249),
  init_LL = fit0$LL,
  init_FF = fit0$FF
)
readr::write_rds(fit_approx, "~/Documents/data/passPCA/lsa_k13_cheby_approx.rds")

init_W <- fit0$LL
init_H <- t(fit0$FF)

fit <- nnmf(
  A = as.matrix(Y_tilde),
  init = list(W = init_W, H = init_H),
  k = 13
)

readr::write_rds(fit, "~/Documents/data/passPCA/lsa_k13_frob_fit.rds")

library(fastTopics)
structure_plot(
  log1pNMF:::normalize_bars(fit$W)[i,],
  grouping = clusters[i],
  gap = 20
)

colnames(fit$W) <- NULL
structure_plot(
  log1pNMF:::normalize_bars(fit$W)[i,],
  grouping = conditions[i],
  gap = 20
)

celltype_factors <- c(
  1, 2, 3, 5, 6, 7, 8, 9, 12
)
other_factors <- c(
  4, 10, 11, 12, 13
)

structure_plot(
  log1pNMF:::normalize_bars(fit$W)[i, celltype_factors],
  grouping = clusters[i],
  gap = 20
)

structure_plot(
  log1pNMF:::normalize_bars(fit$W)[i, other_factors],
  grouping = conditions[i],
  gap = 20
)

load("../data/experiment_results.Rdata")

log1p_k13 <- res_list$pancreas$`1`

hoyer_sparsity <- function(x) {
  
  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))
  
}

structure_plot(
  log1pNMF:::normalize_bars(fit_approx$LL)[i,],
  grouping = clusters[i],
  gap = 20
)

celltype_factors_cheb <- c(
  1, 2, 3, 4, 5, 6, 8, 9, 11, 12
)
other_factors_cheb <- c(
  7, 10, 13
)

structure_plot(
  log1pNMF:::normalize_bars(fit_approx$LL)[i, celltype_factors_cheb],
  grouping = clusters[i],
  gap = 20
)

structure_plot(
  log1pNMF:::normalize_bars(fit_approx$LL)[i, other_factors_cheb],
  grouping = conditions[i],
  gap = 20
)

mean(apply(log1p_k13$LL, 2, hoyer_sparsity))
mean(apply(fit$W, 2, hoyer_sparsity))
mean(apply(fit_approx$LL, 2, hoyer_sparsity))
mean(apply(log1p_k13$FF, 2, hoyer_sparsity))
mean(apply(fit_approx$FF, 2, hoyer_sparsity))
mean(apply(t(fit$H), 2, hoyer_sparsity))

