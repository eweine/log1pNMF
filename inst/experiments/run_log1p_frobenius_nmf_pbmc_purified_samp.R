args <- commandArgs(trailingOnly = TRUE)
cc <- as.numeric(args[1])

load("~/Documents/data/fastglmpca/raw_data/pbmc_purified.RData")
set.seed(1)
library(Seurat)

so <- Seurat::CreateSeuratObject(counts = Matrix::t(counts))

so <- FindVariableFeatures(so, nfeatures = ceiling(.25 * ncol(counts)))

idx <- sample(1:nrow(counts), ceiling(.1 * nrow(counts)))

counts_samp <- counts[idx, ]

counts_samp <- counts_samp[, VariableFeatures(so)]

counts_samp <- Matrix::rowScale(counts_samp, 1 / Matrix::rowSums(counts_samp))
counts_samp <- MatrixExtra::mapSparse(counts_samp, log1p)

library(RcppML)

nmf_out <- nmf(counts_samp, k = 10)
readr::write_rds(
  nmf_out,
  "~/Documents/passPCA/inst/experiments/results/log1p_frobenius_nmf_pbmc_purified_samp.rds"
)
