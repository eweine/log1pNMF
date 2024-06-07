# Here, I want to input my code and try and reproduce the plots that Jason made
# I will start with just one
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
log1p_fit <- readr::read_rds(
  "~/Documents/passPCA/inst/experiments/results/log1p_quad_approx_pbmc_purified_30_factors_100_iter_nmf_init.rds"
)
pois_nmf_fit <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/pois_nmf_pbmc_purified_30_factors_125_iter.rds")
frob_nmf_fit <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/log1p_transformation_pbmc_purified_frob_nmf.rds")

load("~/Documents/data/fastglmpca/raw_data/pbmc_purified.RData")

FF_log1p <- log1p_fit$U
FF_log1p <- scale(FF_log1p, center = FALSE, scale = apply(FF_log1p, 2, max))
colnames(FF_log1p) <- paste0("k", 1:ncol(FF_log1p))

FF_frob_nmf <- frob_nmf_fit$w * sqrt(frob_nmf_fit$d)
FF_frob_nmf <- scale(FF_frob_nmf, center = FALSE, scale = apply(FF_frob_nmf, 2, max))
colnames(FF_frob_nmf) <- paste0("k", 1:ncol(FF_frob_nmf))

FF_pois_nmf <- pois_nmf_fit$L

cell.type <- samples$celltype

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
FF_frob_nmf <- FF_frob_nmf[cell.idx, ]
FF_pois_nmf <- FF_pois_nmf[cell.idx, ]

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
frob_nmf_tib <- make.heatmap.tib(FF_frob_nmf)
pois_nmf_tib <- make.heatmap.tib(FF_pois_nmf)

heatmap.tib <- log1p_tib %>% mutate(Method = "log1p Poisson NMF") %>%
  bind_rows(frob_nmf_tib %>% mutate(Method = "log1p Transformation Frob. NMF")) %>%
  bind_rows(pois_nmf_tib %>% mutate(Method = "Poisson NMF")) %>%
  mutate(Method = factor(Method, levels = c("log1p Poisson NMF", "log1p Transformation Frob. NMF", "Poisson NMF")))

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
  theme_minimal() +
  geom_hline(yintercept = -cell_type_breaks, size = 0.1) +
  facet_wrap(~Method, ncol = 1) +
  theme(legend.position = "none", 
        strip.text = element_text(size = 16)) 

# Now that I have this, I think it would be good to share it in a workflowr


plt