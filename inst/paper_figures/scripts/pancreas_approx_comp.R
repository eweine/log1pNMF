library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)
library(passPCA)

set.seed(1)

load("~/Downloads/pancreas.RData")

i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]
sample_info <- transform(sample_info,celltype = factor(celltype))

s <- Matrix::rowSums(counts)
s <- s / mean(s)

genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]

n <- nrow(counts)
p <- ncol(counts)

K <- 9

# now, I want to fit (1) the approximate loglik approach and (2) frob nmf
set.seed(1)
fit_exact_5iter <- fit_poisson_log1p_nmf(
  Y = counts,
  K = K,
  s = s,
  cc = 1,
  loglik = "exact",
  init_method = "rank1",
  control = list(
    maxiter = 5,
    threads = 7
  )
)

fit_exact <- fit_poisson_log1p_nmf(
  Y = counts,
  init_LL = fit_exact_5iter$LL,
  init_FF = fit_exact_5iter$FF,
  s = s,
  cc = 1,
  loglik = "exact",
  control = list(
    maxiter = 100,
    threads = 7
  )
)

fit_approx <- fit_poisson_log1p_nmf(
  Y = counts,
  init_LL = fit_exact_5iter$LL,
  init_FF = fit_exact_5iter$FF,
  s = s,
  cc = 1,
  loglik = "approx",
  control = list(
    maxiter = 100,
    threads = 7
  )
)

counts <- Matrix::Diagonal(x = 1 / s) %*% counts
counts <- MatrixExtra::mapSparse(counts, log1p)

library(NNLM)
set.seed(1)

init_list <- list(
  W = fit_exact_5iter$LL,
  H = t(fit_exact_5iter$FF)
)

nmf_k9 <- nnmf(
  A = as.matrix(counts),
  init = init_list,
  k = 9
)

celltype <- sample_info$celltype
celltype <-
  factor(celltype,
         c("acinar","ductal","activated_stellate","quiescent_stellate",
           "endothelial","macrophage","mast","schwann","alpha","beta",
           "delta","gamma","epsilon"))

levels(celltype)[levels(celltype) == "activated_stellate"] <- "activated PSC"
levels(celltype)[levels(celltype) == "quiescent_stellate"] <- "quiescent PSC"

set.seed(1)
colnames(fit_exact$LL) <- paste0("k", 1:9)
colnames(fit_approx$LL) <- paste0("k", 1:9)

topic_order <- paste0("k",
                      c(1, 5, 7, 8, 3, 4, 9, 6, 2))

log1p_sp <- normalized_structure_plot(
  fit_exact,
  topics = topic_order,
  grouping = celltype,gap = 30,perplexity = 70,n = Inf, font.size = 12
) 

g_log1p_sp <- log1p_sp$plot + 
  ggtitle("Loadings of log1p Link Poisson NMF With c = 1 - Exact Optimization") +
  theme(
    plot.title = element_text(size = 10)
  ) + ylim(0, 2.6)

g_log1p_sp_approx <- normalized_structure_plot(
  fit_approx, 
  topics = topic_order,
  grouping = celltype,
  gap = 30,perplexity = 70, font.size = 12, loadings_order = log1p_sp$loadings_order
)$plot + ggtitle("Loadings of log1p Link Poisson NMF With c = 1 - Approximate Optimization") +
  theme(
    plot.title = element_text(size = 10)
  ) + ylim(0, 2.6)

colnames(nmf_k9$W) <- paste0(
  "k", c(1, 2, 4, 3, 5, 6, 7, 8, 9)
)

nmf_k9$W <- nmf_k9$W[,paste0("k", 1:9)]

g_frob_sp <- structure_plot(
  passPCA:::normalize_bars(nmf_k9$W), 
  grouping = celltype,
  topics = topic_order,
  gap = 30,perplexity = 70, font.size = 12, loadings_order = log1p_sp$loadings_order
)$plot + ggtitle("Loadings of Frobenius NMF Applied to log1p Transformed Data") +
  theme(
    plot.title = element_text(size = 10)
  ) + ylim(0, 2.6)

cor_mat <- cor(fit_exact$FF, method = "spearman")
cor_exact <- median(cor_mat[lower.tri(cor_mat)])

cor_mat2 <- cor(fit_approx$FF, method = "spearman")
cor_approx <- median(cor_mat2[lower.tri(cor_mat2)])

cor_mat3 <- cor(t(nmf_k9$H), method = "spearman")
cor_frob <- median(cor_mat3[lower.tri(cor_mat3)])

hoyer_sparsity <- function(x) {
  
  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))
  
}

lsp_exact <- median(apply(
  fit_exact$LL, 2, hoyer_sparsity
))

fsp_exact <- median(apply(
  fit_exact$FF, 2, hoyer_sparsity
))

lsp_approx <- median(apply(
  fit_approx$LL, 2, hoyer_sparsity
))

fsp_approx <- median(apply(
  fit_approx$FF, 2, hoyer_sparsity
))

lsp_frob <- median(apply(
  nmf_k9$W, 2, hoyer_sparsity
))

fsp_frob <- median(apply(
  t(nmf_k9$H), 2, hoyer_sparsity
))

# 1) Correlation
df_correlation <- data.frame(
  Method      = c("Exact", "Approx", "Frob"),
  Correlation = c(cor_exact, cor_approx, cor_frob)
)

p1 <- ggplot(df_correlation, aes(x = Method, y = Correlation)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Median Factor Correlation") +
  cowplot::theme_cowplot()

# 2) L-sparsity
df_l_sparsity <- data.frame(
  Method     = c("Exact", "Approx", "Frob"),
  L_Sparsity = c(lsp_exact, lsp_approx, lsp_frob)
)

p2 <- ggplot(df_l_sparsity, aes(x = Method, y = L_Sparsity)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Median Loading Sparsity") +
  cowplot::theme_cowplot()

# 3) F-sparsity
df_f_sparsity <- data.frame(
  Method     = c("Exact", "Approx", "Frob"),
  F_Sparsity = c(fsp_exact, fsp_approx, fsp_frob)
)

p3 <- ggplot(df_f_sparsity, aes(x = Method, y = F_Sparsity)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Median Factor Sparsity") +
  cowplot::theme_cowplot()

g1 <- ggarrange(
  g_log1p_sp, g_log1p_sp_approx, g_frob_sp,
  nrow = 3, ncol = 1,
  labels = c("A", "B", "C")
)

g2 <- ggarrange(
  p1, p2, p3,
  nrow = 3, ncol = 1,
  labels = c("D", "E", "F")
)

g <- ggarrange(g1, g2, nrow = 1, ncol = 2, widths = c(2.5, 1))

ggsave(
  "/Users/eweine/Documents/passPCA/inst/paper_figures/pdfs/pancreas_approx_comp.pdf",
  g,
  device = "pdf",
  width = 11,
  height = 9
)

