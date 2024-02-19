args <- commandArgs(trailingOnly = TRUE)
cc <- as.numeric(args[1])
#cc <- 1

load("~/Documents/data/fastglmpca/raw_data/pbmc_purified.RData")
set.seed(1)
library(Seurat)

so <- Seurat::CreateSeuratObject(counts = Matrix::t(counts))

so <- FindVariableFeatures(so, nfeatures = ceiling(.25 * ncol(counts)))

idx <- sample(1:nrow(counts), ceiling(.1 * nrow(counts)))

counts_samp <- counts[idx, ]

counts_samp <- counts_samp[, VariableFeatures(so)]

counts_samp <- counts_samp[, Matrix::colSums(counts_samp) > 0]

rm(counts)
rm(genes)
rm(samples)
gc()

library(passPCA)

log1p_fit <- fit_factor_model_log1p_exact(
  Y = counts_samp,
  K = 10,
  maxiter = 1000,
  s = cc * as.vector(Matrix::rowSums(counts_samp) / mean(Matrix::rowSums(counts_samp)))
)

readr::write_rds(
  log1p_fit,
  glue::glue(
    "results/log1p_exact_c{cc}_pbmc_purified_10_factors_1000_iter.rds"
  )
)
