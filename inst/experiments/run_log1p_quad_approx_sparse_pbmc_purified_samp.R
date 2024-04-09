load("~/Documents/data/fastglmpca/raw_data/pbmc_purified.RData")
set.seed(1)
counts <- counts[, Matrix::colSums(counts) > 0]
sz <- Matrix::rowSums(counts) / mean(Matrix::rowSums(counts))

library(Seurat)
colnames(counts) <- gsub("_", "-", genes$symbol)
nd <- colnames(counts)
nd <- nd[!duplicated(nd)]
counts <- counts[, nd]

so <- CreateSeuratObject(Matrix::t(counts))

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <= 5)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 4000)

counts <- counts[rownames(so@meta.data), ]
sz <- sz[rownames(so@meta.data)]
counts <- counts[, VariableFeatures(so)]
counts <- counts[, Matrix::colSums(counts) >= 100]


library(passPCA)

set.seed(1)
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 30,
  maxiter = 50,
  approx_range = c(0, 1.25),
  s = as.vector(sz),
  init_method = "frob_nmf"
)

rownames(log1p_fit$V) <- colnames(counts)
rownames(log1p_fit$U) <- rownames(counts)


library(Matrix)
set.seed(1)
log1p_fit_frob_init <- fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 30,
  maxiter = 100,
  approx_range = c(0, 1.25),
  s = as.vector(Matrix::rowSums(counts) / mean(Matrix::rowSums(counts))),
  init_method = "frob_nmf"
)

U_normalized <- apply(log1p_fit_frob_init$U, 2, function(col) col / max(col))

fastTopics::structure_plot(
  U_normalized,
  grouping = samples$celltype,
  n = 5000
)

mod2 <- readr::read_rds(
  glue::glue(
    "results/log1p_quad_approx_pbmc_purified_30_factors_125_iter.rds"
  )
)

U_normalized <- t(apply(mod2$U, 1, function(row) row / sum(row)))


readr::write_rds(
  log1p_fit,
  glue::glue(
    "results/log1p_quad_approx_pbmc_purified_30_factors_125_iter.rds"
    )
)

library(fastTopics)
set.seed(1)
nmf_fit <- fit_poisson_nmf(
  X = counts,
  k = 30,
  numiter = 125
)

readr::write_rds(
  nmf_fit,
  glue::glue(
    "results/pois_nmf_pbmc_purified_30_factors_125_iter.rds"
  )
)

sz <- as.vector(Matrix::rowSums(counts) / mean(Matrix::rowSums(counts)))

counts <- counts / sz
counts <- MatrixExtra::mapSparse(counts, log1p)

frob_nmf <- RcppML::nmf(counts, k = 30)

U_nmf <- frob_nmf$w * sqrt(frob_nmf$d)
V_nmf <- sqrt(frob_nmf$d) * frob_nmf$h

U_normalized_nmf <- apply(U_nmf, 1, function(row) row / sum(row))

readr::write_rds(frob_nmf, "results/log1p_transformation_pbmc_purified_frob_nmf.rds")
