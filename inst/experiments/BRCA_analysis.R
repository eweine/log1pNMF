library(Matrix)
library(tidyr)
library(stringr)
library(purrr)

dat1 <- readr::read_rds("~/Downloads/1863-counts_cells_cohort1.rds")
md1 <- readr::read_csv("~/Downloads/1872-BIOKEY_metaData_cohort1_web.csv")

dat2 <- readr::read_rds("~/Downloads/1867-counts_cells_cohort2.rds")
md2 <- readr::read_csv("~/Downloads/1871-BIOKEY_metaData_cohort2_web.csv")

common_genes <- intersect(rownames(dat1), rownames(dat2))

dat1 <- Matrix::t(dat1[common_genes, ])
dat2 <- Matrix::t(dat2[common_genes, ])
gc()

library(dplyr)
md1 <- md1 %>% dplyr::filter(
  cellType %in% c("T_cell", "B_cell", "Myeloid_cell")
)
md2 <- md2 %>% dplyr::filter(
  cellType %in% c("T_cell", "B_cell", "Myeloid_cell")
)

dat1 <- dat1[md1$Cell, ]
dat2 <- dat2[md2$Cell, ]

high_genes <- names(which(Matrix::colSums(dat1) + Matrix::colSums(dat2) > 0))
dat1 <- dat1[, high_genes]
dat2 <- dat2[, high_genes]

nrd1 <- nrow(dat1)

m1 <- Matrix::summary(dat1)
m2 <- Matrix::summary(dat2)

dn_list <- list()
dn_list[[1]] <- c(rownames(dat1), rownames(dat2))
dn_list[[2]] <- colnames(dat1)

rm(dat1, dat2, md1, md2, common_genes, high_genes)
gc()

counts <- sparseMatrix(
  i = c(
    m1$i,
    m2$i + nrd1
  ),
  j = c(
    m1$j,
    m2$j
  ),
  x = c(
    m1$x,
    m2$x
  ),
  dimnames = dn_list
)

rm(dn_list, m1, m2, nrd1)
gc()

library(Seurat)

rs <- Matrix::rowSums(counts)
s <- rs / mean(rs)

counts <- Matrix::t(counts)
rm(rs)
gc()

counts <- counts[!grepl("\\.", rownames(counts)), ]

so <- CreateSeuratObject(counts = counts)

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

med_nFeature <- median(so[["nFeature_RNA"]][,1])
mad_n_Feature <- mad(so[["nFeature_RNA"]][,1])

so <- subset(so, subset = (nFeature_RNA > 200) & (nFeature_RNA < med_nFeature + 2 * mad_n_Feature) & percent.mt <= 10)

so <- NormalizeData(so)

so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 500)

counts <- Matrix::t(counts)
counts <- counts[colnames(so@assays$RNA$counts), ]
counts <- counts[, VariableFeatures(so)]
counts <- counts[, Matrix::colSums(counts) >= 100]

rm(so, mad_n_Feature, med_nFeature)
gc()

counts <- counts[Matrix::rowSums(counts) > 0, ]

set.seed(1)
log1p_fit <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  counts,
  K = 20,
  approx_range = c(0, 1.25),
  maxiter = 50,
  init_method = "frob_nmf",
  s = s[rownames(counts)]
)

rownames(log1p_fit$V) <- colnames(counts)
rownames(log1p_fit$U) <- rownames(counts)

#readr::write_rds(log1p_fit, "~/Documents/passPCA/inst/experiments/results/brca_immune_cell_log1p.rds")
log1p_fit <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/brca_immune_cell_log1p.rds")

#rm(log1p_fit)
#gc()

# set.seed(1)
# log1p_fit_greedy <- passPCA::fit_factor_model_log1p_quad_approx_sparse_greedy(
#   counts,
#   K = 20,
#   approx_range = c(0, 1.25),
#   iter_per_factor = 10,
#   s = s[rownames(counts)]
# )

#rownames(log1p_fit_greedy$V) <- colnames(counts)
#rownames(log1p_fit_greedy$U) <- rownames(counts)

#readr::write_rds(log1p_fit_greedy, "~/Documents/passPCA/inst/experiments/results/brca_immune_cell_log1p_greedy.rds")
log1p_fit_greedy <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/brca_immune_cell_log1p_greedy.rds")

set.seed(1)
# log1p_fit_greedy_init <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
#   counts,
#   K = 20,
#   approx_range = c(0, 1.25),
#   maxiter = 50,
#   init_U = log1p_fit_greedy$U,
#   init_V = log1p_fit_greedy$V,
#   s = s[rownames(counts)]
# )

#rownames(log1p_fit_greedy_init$V) <- colnames(counts)
#rownames(log1p_fit_greedy_init$U) <- rownames(counts)

#readr::write_rds(log1p_fit_greedy_init, "~/Documents/passPCA/inst/experiments/results/brca_immune_cell_log1p_greedy_init.rds")
log1p_fit_greedy_init <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/brca_immune_cell_log1p_greedy_init.rds")

md1 <- readr::read_csv("~/Downloads/1872-BIOKEY_metaData_cohort1_web.csv")
md2 <- readr::read_csv("~/Downloads/1871-BIOKEY_metaData_cohort2_web.csv")

md <- rbind(md1, md2)

md <- md %>%
  dplyr::filter(
    Cell %in% rownames(counts)
  )

md$ct_treat <- paste(md$cohort, md$cellType, md$timepoint, sep = "-")

s_df <- md %>% dplyr::group_by(ct_treat) %>%
  summarize(
    tot = n()
  )

# Now, I want to make a heatmap...
FF_log1p <- log1p_fit$U
FF_log1p <- scale(FF_log1p, center = FALSE, scale = apply(FF_log1p, 2, max))
colnames(FF_log1p) <- paste0("k", 1:ncol(FF_log1p))

# FF_log1p_greedy <- log1p_fit_greedy$U
# FF_log1p_greedy <- scale(
#   FF_log1p_greedy,
#   center = FALSE,
#   scale = apply(FF_log1p_greedy, 2, max)
# )
# colnames(FF_log1p_greedy) <- paste0("k", 1:ncol(FF_log1p_greedy))
#
# FF_log1p_greedy_init <- log1p_fit_greedy_init$U
# FF_log1p_greedy_init <- scale(
#   FF_log1p_greedy_init,
#   center = FALSE,
#   scale = apply(FF_log1p_greedy_init, 2, max)
# )
# colnames(FF_log1p_greedy_init) <- paste0("k", 1:ncol(FF_log1p_greedy_init))


cell.type <- factor(md$ct_treat)

# Downsample the number of cells and sort them using tSNE.
set.seed(8675309)
cell.idx <- numeric(0)
cell.types <- levels(cell.type)
for (i in 1:length(cell.types)) {
  which.idx <- which(cell.type == cell.types[i])
  # Downsample common cell types.
  if (length(which.idx) > 1000) {
    which.idx <- sample(which.idx, 1000)
  }
  # Don't include rare cell types.
  if (length(which.idx) > 20) {
    # Sort using tsne.
    tsne.res <- Rtsne::Rtsne(
      FF_log1p[which.idx, ],
      dims = 1,
      pca = FALSE,
      normalize = FALSE,
      perplexity = min(100, floor((length(which.idx) - 1) / 3) - 1),
      theta = 0.1,
      max_iter = 1000,
      eta = 200,
      check_duplicates = FALSE
    )$Y[, 1]
    which.idx <- which.idx[order(tsne.res)]
    cell.idx <- c(cell.idx, which.idx)
  }
}

cell.type <- cell.type[cell.idx]
cell.type <- droplevels(cell.type)

FF_log1p <- FF_log1p[cell.idx, ]
#FF_log1p_greedy <- FF_log1p_greedy[cell.idx, ]
#FF_log1p_greedy_init <- FF_log1p_greedy_init[cell.idx, ]

make.heatmap.tib <- function(FF) {
  tib <- as_tibble(scale(FF, center = FALSE, scale = apply(FF, 2, max))) %>%
    mutate(Cell.type = cell.type) %>%
    arrange(Cell.type) %>%
    mutate(Cell.idx = row_number())

  tib <- tib %>%
    pivot_longer(
      -c(Cell.idx, Cell.type),
      names_to = "Factor",
      values_to = "Loading",
      values_drop_na = TRUE
    ) %>%
    mutate(Factor = as.numeric(str_extract(Factor, "[0-9]+")))

  return(tib)
}

log1p_tib <- make.heatmap.tib(FF_log1p)
#log1p_greedy_tib <- make.heatmap.tib(FF_log1p_greedy)
#log1p_greedy_init_tib <- make.heatmap.tib(FF_log1p_greedy_init)

heatmap.tib <- log1p_tib %>% mutate(Method = "log1p Poisson NMF") %>%
  mutate(Method = factor(Method, levels = c("log1p Poisson NMF")))
  #bind_rows(log1p_greedy_tib %>% mutate(Method = "log1p Poisson NMF greedy")) %>%
  #bind_rows(log1p_greedy_init_tib %>% mutate(Method = "log1p Poisson NMF greedy + backfit")) %>%
  #mutate(Method = factor(Method, levels = c("log1p Poisson NMF", "log1p Poisson NMF with rank-1 initialization", "log1p Transformation Frob. NMF", "Poisson NMF")))

tib <- heatmap.tib %>%
  group_by(Cell.type, Cell.idx) %>%
  summarize()

cell_type_breaks <- c(1, which(tib$Cell.type[-1] != tib$Cell.type[-nrow(tib)]))
label_pos <- cell_type_breaks / 2 + c(cell_type_breaks[-1], nrow(tib)) / 2

library(ggplot2)

plt <- ggplot(heatmap.tib, aes(x = Factor, y = -Cell.idx, fill = Loading)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(y = "") +
  scale_y_continuous(breaks = -label_pos,
                     minor_breaks = NULL,
                     labels = levels(cell.type)) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  theme_minimal() +
  geom_hline(yintercept = -cell_type_breaks, size = 0.1) +
  facet_wrap(~Method, ncol = 1, axes = "all") +
  theme(legend.position = "none",
        strip.text = element_text(size = 16))

plt

