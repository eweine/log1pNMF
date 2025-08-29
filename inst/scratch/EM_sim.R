library(dplyr)
library(fastTopics)
library(log1pNMF)
library(ggplot2)
library(cowplot)
library(Matrix)
library(readr)

set.seed(1)

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

i <- which(clusters == "Endothelial/Mesnchymal")

tm_k13 <- res_list$pancreas$`Inf`

#L <- tm_k13$L
colnames(tm_k13$L) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
tm_k13$L <- tm_k13$L[i,c("k9", "k12")]

colnames(tm_k13$F) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
tm_k13$F <- tm_k13$F[,c("k9", "k12")]

L <- poisson2multinom(tm_k13)$L
FF <- poisson2multinom(tm_k13)$F

Pi <- L %*% t(FF)

s <- Matrix::rowSums(counts[i,])

Y_sim <- matrix(
  data = 0,
  nrow = length(s),
  ncol = ncol(Pi)
)

for (j in 1:length(s)) {
  
  Y_sim[j, ] <- rmultinom(n = 1, size = s[j], prob = Pi[j,])
  
}


genes_to_use <- which(Matrix::colSums(Y_sim>0)>=2)
Y_sim <- Y_sim[, genes_to_use]

Y_sim <- as(Y_sim, "CsparseMatrix")

library(fastTopics)

set.seed(1)
ft_fit <- fit_poisson_nmf(X = Y_sim, k = 2)

structure_plot(ft_fit)

library(log1pNMF)

set.seed(1)
log1p_fit <- fit_poisson_log1p_nmf(
  Y = Y_sim, K = 2, loglik = "exact", init_method = "random"
)

normalized_structure_plot(log1p_fit)
