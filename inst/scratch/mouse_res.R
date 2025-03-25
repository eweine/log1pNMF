library(dplyr)
load("~/Downloads/mouse_brain_stim.Rdata")

cells <- cells %>% dplyr::filter(!is.na(maintype))

counts <- counts[rownames(counts) %in% cells$`...1`, ]
counts <- counts[, Matrix::colSums(counts) > 0]
counts <- counts[Matrix::rowSums(counts) > 0, ]

#nmf <- readr::read_rds("/Users/eweine/Documents/data/mouse_log1p_c1e-04_k15_exact_100_iter.rds")

library(passPCA)

#normalized_structure_plot(mod, grouping = cells$maintype)

library(SingleCellExperiment)
cells <- as.data.frame(cells)
rownames(cells) <- cells$`...1`
cells$sample <- sub("^(([^_]*)_([^_]*)).*", "\\1", cells$sample)
cells <- cells %>% dplyr::select(-`...1`)

sce <- SingleCellExperiment(
  assays = list(counts = Matrix::t(counts)),  # Add the counts matrix
  colData = cells                    # Add metadata
)

sb <- pseudobulk(sce, group_by = vars(maintype, sample, stim))
sb_f <- sb[,sb$maintype == "Excitatory"]


mod <- glm_gp(data = sb_f, design = ~1 + stim)
de <- test_de(mod, contrast = "stim4h")



# now, I want to look at the actual effect sizes
# but first I want to do some filtering

b <- data.frame(mod$Beta)
b$gene <- rownames(b)

b <- b %>%
  dplyr::filter(Intercept > -1e5)

de <- de %>%
  dplyr::filter(
    abs(lfc) < 10
  )

b <- b %>%
  dplyr::filter(rownames(b) %in% de$name)

b <- b %>%
  dplyr::filter(stim1h > -1e5)

b <- b %>%
  dplyr::mutate(
    ctrl_expr = exp(Intercept),
    stim1_expr = exp(Intercept + stim1h),
    stim4_expr = exp(Intercept + stim4h)
  )

b <- b %>%
  dplyr::mutate(
    stim1_diff = stim1_expr - ctrl_expr,
    stim4_diff = stim4_expr - ctrl_expr,
    stim1_log1p_diff = log1p(stim1_expr) - log1p(ctrl_expr),
    stim4_log1p_diff = log1p(stim4_expr) - log1p(ctrl_expr)
  )

de_high <- de %>% dplyr::filter(adj_pval < 0.1)

b_pos <- b %>%
  dplyr::filter(rownames(b) %in% de_high$name) %>%
  dplyr::filter(stim1_diff > 0) 


