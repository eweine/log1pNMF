library(Seurat)
library(magrittr)
so <- SeuratDisk::LoadH5Seurat("~/Downloads/GSE226824_HSPC-all_filtered.h5seurat")

counts <- Matrix::t(so@assays$RNA@counts)
counts <- counts[,Matrix::colSums(counts) > 0]
gc()

rs <- Matrix::rowSums(counts)
s <- rs / mean(rs)

set.seed(1)
log1p_k1 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 1,
  approx_range = c(0, 1.25),
  maxiter = 10,
  s = s
)

n <- nrow(counts)
p <- ncol(counts)
K <- 11

init_LL <- log1p_k1$U %>%
  cbind(
    matrix(
      data = rexp(
        n = n * (K - 1), rate = 15
      ),
      nrow = n,
      ncol = K - 1
    )
  )

init_FF <- log1p_k1$V %>%
  cbind(
    matrix(
      data = rexp(
        n = p * (K - 1), rate = 15
      ),
      nrow = p,
      ncol = K - 1
    )
  )

log1p_k1 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = K,
  maxiter = 100,
  approx_range = c(0, 1.25),
  init_U = init_LL,
  init_V = init_FF,
  s = s
)

normalize_bars <- function(LL) {

  max_col <- apply(LL, 2, max)
  sweep(LL, 2, max_col, FUN = "/")

}

readr::write_rds(log1p_k1, "~/Documents/data/passPCA/experiment_results/log1p_k11_cytokine_stim_ts.rds")

LL <- log1p_k1$U
LL <- normalize_bars(LL)

library(fastTopics)
structure_plot(LL[,2:11], grouping = paste0(sub("#.*", "", so@meta.data$clusters), as.numeric(so@meta.data$time)), gap = 20)


paste0(sub("#.*", "", so@meta.data$clusters), so@meta.data$time)
library(RcppML)

# would also be interesting to look at cNMF fit and other fits
