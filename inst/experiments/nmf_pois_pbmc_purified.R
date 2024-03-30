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

#rm(counts)
#rm(genes)
samples <- samples[idx, ]
gc()

library(fastTopics)

#nmf_pois <- fit_poisson_nmf(
#  counts_samp,
#  k = 10,
#  numiter = 1000
#)

#readr::write_rds(
#  nmf_pois,
#  glue::glue(
#    "results/pois_nmf_pbmc_purified_10_factors_1000_iter.rds"
#  )
#)

nmf_pois <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/pois_nmf_pbmc_purified_10_factors_1000_iter.rds")

library(fastTopics)

structure_plot(
  fit = nmf_pois,
  grouping = samples$celltype,
  gap = 20
)
library(fastTopics)
log1p_fits <- readr::read_rds(
  "~/Documents/passPCA/inst/experiments/results/log1p_exact_c0.1_pbmc_purified_10_factors_1000_iter.rds"
)

log1p_fits$U <- t(t(log1p_fits$U) / apply(log1p_fits$U,2,max))


structure_plot(fit = log1p_fits$U,
               grouping = samples$celltype,
               gap = 20, n = 5000)
