library(Matrix)
library(log1pNMF)

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")

library(NNLM)

K <- 13

s <- Matrix::rowSums(counts) 
s <- s / mean(s)
Y_tilde <- Matrix::Diagonal(x = 1/s) %*% counts
Y_tilde <- MatrixExtra::mapSparse(Y_tilde, log1p)

set.seed(1)
fit0 <- nnmf(
  A = as.matrix(Y_tilde)
)

init_W <- cbind(
  fit0$W,
  matrix(
    data = 1e-5,
    nrow = nrow(counts),
    ncol = K - 1
  )
)


init_H <- rbind(
  fit0$H,
  matrix(
    data = 1e-5,
    nrow = K - 1,
    ncol = ncol(counts)
  )
)

fit <- nnmf(
  A = as.matrix(Y_tilde),
  init = list(W = init_W, H = init_H),
  k = 13
)

clusters   <- factor(barcodes$celltype,
                     c("Acinar","Ductal","Endothelial/Mesnchymal","Macrophage",
                       "Alpha","Beta","Delta","Gamma"))

library(fastTopics)
structure_plot(
  log1pNMF:::normalize_bars(fit$W),
  grouping = clusters,
  gap = 20
)
