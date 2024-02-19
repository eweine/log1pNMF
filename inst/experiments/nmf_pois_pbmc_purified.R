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

library(fastTopics)

nmf_pois <- fit_poisson_nmf(
  counts_samp,
  k = 10,
  numiter = 1000
)

readr::write_rds(
  nmf_pois,
  glue::glue(
    "results/pois_nmf_pbmc_purified_10_factors_1000_iter.rds"
  )
)


