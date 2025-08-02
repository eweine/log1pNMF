load("~/Documents/data/passPCA/liu_data.Rdata")

genes_to_use <- which(Matrix::colSums(counts>0)>=3)

counts <- counts[,genes_to_use]

library(log1pNMF)
set.seed(1)
log1p_exact_mod <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 15,
  loglik = "exact"
)

library(log1pNMF)
set.seed(1)
log1p_approx_mod <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 15,
  loglik = "approx",
  approx_technique = "chebyshev"
)

# normalized_structure_plot(
#   log1p_exact_mod,
#   grouping = ct$celltype,
#   gap = 10
# )

s <- Matrix::rowSums(counts)
s <- s / mean(s)

Y_tilde <- Matrix::Diagonal(x = 1/s) %*% counts

Y_tilde <- MatrixExtra::mapSparse(Y_tilde, log1p)

Y_dense <- as.matrix(Y_tilde)

set.seed(1)
nnlm_r1 <- NNLM::nnmf(
  A = Y_dense,
  k = 1
)

L_nnlm_init <- cbind(
  nnlm_r1$W,
  matrix(
    data = 1e-5,
    nrow = nrow(nnlm_r1$W),
    ncol = 14
  )
)

F_nnlm_init <- rbind(
  nnlm_r1$H,
  matrix(
    data = 1e-5,
    ncol = ncol(nnlm_r1$H),
    nrow = 14
  )
)

set.seed(1)
nnlm_full <- NNLM::nnmf(
  A = Y_dense,
  k = 15,
  init = list(W = L_nnlm_init, H = F_nnlm_init)
)

structure_plot(
  log1pNMF:::normalize_bars(nnlm_full$W), 
  grouping = ct$celltype,gap = 25,perplexity = 70, font.size = 14
)

normalized_structure_plot(
  log1p_approx_mod,
  grouping = ct$celltype,
  gap = 25
)

hoyer_sparsity <- function(x) {
  
  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))
  
}
