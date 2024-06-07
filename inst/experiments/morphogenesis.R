library(dplyr)
counts <- Matrix::readMM(
  "~/Downloads/GSE175634_cell_counts.mtx"
)

cell_md <- readr::read_tsv("~/Downloads/GSE175634_cell_metadata.tsv")

s_df <- cell_md %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(
    tot = n()
  )

cell_md <- cell_md %>%
  dplyr::filter(
    individual == 19093
  )

cell_idx <- readr::read_tsv("~/Downloads/GSE175634_cell_indices.tsv")
gene_idx <- readr::read_tsv("~/Downloads/GSE175634_gene_indices_counts.tsv")

rownames(counts) <- gene_idx$gene_name
colnames(counts) <- cell_idx$cell_name
counts <- Matrix::t(counts)

cell_md <- cell_md %>% dplyr::filter(
  type != "UNK"
)

counts <- counts[cell_md$cell, ]
counts <- counts[,Matrix::colSums(counts) > 10]

rm(cell_idx, gene_idx, s_df)
gc()

library(passPCA)

rs <- Matrix::rowSums(counts)
s <- rs / mean(rs)

rm(rs)
gc()

counts <- as(counts, "CsparseMatrix")

library(Matrix)
set.seed(1)
mod2 <- fit_factor_model_log1p_quad_approx_sparse_greedy(
  Y = counts,
  s = s,
  K = 30,
  iter_per_factor = 10,
  approx_range = c(0, 1.25)
)

readr::write_rds(mod, "~/Documents/data/passPCA/experiment_results/log1p_morphogenesis_K50.rds")

cell_md$diffday <- factor(cell_md$diffday, levels = c("day0", "day1", "day3", "day5", "day7", "day11", "day15"))

FF_log1p_greedy <- mod$U
FF_log1p_greedy <- scale(
  FF_log1p_greedy,
  center = FALSE,
  scale = apply(FF_log1p_greedy, 2, max)
)
colnames(FF_log1p_greedy) <- paste0("k", 1:ncol(FF_log1p_greedy))

cell.type <- cell_md$diffday

# Downsample the number of cells and sort them using tSNE.
set.seed(8675309)
cell.idx <- numeric(0)
cell.types <- levels(cell.type)
for (i in 1:length(cell.types)) {
  which.idx <- which(cell.type == cell.types[i])
  # Downsample common cell types.
  if (length(which.idx) > 1250) {
    which.idx <- sample(which.idx, 1250)
  }
  # Don't include rare cell types.
  if (length(which.idx) > 20) {
    # Sort using tsne.
    tsne.res <- Rtsne::Rtsne(
      FF_log1p_greedy[which.idx, ],
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

FF_log1p_greedy <- FF_log1p_greedy[cell.idx, ]

library(plyr)
library(purrr)
library(dplyr)
library(tidyr)
library(stringi)
library(stringr)
make.heatmap.tib <- function(FF) {
  tib <- as_tibble(scale(FF, center = FALSE, scale = apply(FF, 2, max))) %>%
    dplyr::mutate(Cell.type = cell.type) %>%
    dplyr::arrange(Cell.type) %>%
    dplyr::mutate(Cell.idx = dplyr::row_number())

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

log1p_greedy_tib <- make.heatmap.tib(FF_log1p_greedy)

heatmap.tib <- log1p_greedy_tib %>% mutate(Method = "log1p Poisson NMF Greedy Fit") %>%
  mutate(Method = factor(Method, levels = c("log1p Poisson NMF Greedy Fit")))

tib <- heatmap.tib %>%
  dplyr::group_by(Cell.type, Cell.idx) %>%
  dplyr::summarize()



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

