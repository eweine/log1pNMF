load("~/Downloads/mouse_brain_stim.Rdata")

library(magrittr)
library(Matrix)
cells <- cells %>% dplyr::filter(!is.na(maintype))

counts <- counts[cells$...1, ]

counts <- counts[, Matrix::colSums(counts) > 0]

library(passPCA)

s <- Matrix::rowSums(counts)
s <- s/mean(s)

gc()

set.seed(1)
#log1p_mod <- fit_factor_model_log1p_quad_approx_sparse(
#  Y = counts, K = 30, s = s, approx_range = c(0, 1.25), maxiter = 100
#)

counts <- counts / s
counts <- MatrixExtra::mapSparse(counts, log1p)

frob_nmf_mod <- RcppML::nmf(
  A = counts, k = 30
)

d <- sqrt(frob_nmf_mod$d)
LL <- frob_nmf_mod$w * d
FF <- t(d * frob_nmf_mod$h)

#log1p_mod <- readr::read_rds("~/Documents/data/passPCA/experiment_results/mouse_k30.rds")

rownames(LL) <- rownames(counts)
rownames(FF) <- colnames(counts)

# Now, it would be nice to make a heatmap plot
# and output the driving genes for each factor

# It would probably also be useful if I considered different
# initializations
max_col <- apply(LL, 2, max)
log1p_LL <- sweep(LL, 2, max_col, FUN = "/")
log1p_FF <- sweep(FF, 2, max_col, FUN = "*")

#readr::write_rds(log1p_mod, "~/Documents/data/passPCA/experiment_results/mouse_k30.rds")


cells$ct_stim <- factor(paste(cells$maintype, cells$stim))





cell.type <- cells$ct_stim

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
      log1p_LL[which.idx, ],
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


log1p_LL <- log1p_LL[cell.idx, ]

library(tidyr)
library(purrr)
library(dplyr)
library(stringi)
library(stringr)

make.heatmap.tib <- function(LL) {
  tib <- as_tibble(scale(LL, center = FALSE, scale = apply(LL, 2, max))) %>%
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

log1p_tib <- make.heatmap.tib(log1p_LL)
#log1p_k1_init_tib <- make.heatmap.tib(FF_log1p_K1_init)
#frob_nmf_tib <- make.heatmap.tib(FF_frob_nmf)
#pois_nmf_tib <- make.heatmap.tib(FF_pois_nmf)

heatmap.tib <- log1p_tib %>% mutate(Method = "log1p Poisson NMF") %>%
  mutate(Method = factor(Method, levels = c("log1p Poisson NMF")))

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
