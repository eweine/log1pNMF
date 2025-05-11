library(dplyr)
library(ggplot2)
library(ggrepel)
load("../data/raw_data/pancreas.RData")

i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]
sample_info <- transform(sample_info,celltype = factor(celltype))

counts <- counts[, Matrix::colSums(counts) > 0]
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]

load("../data/experiment_results.Rdata")

acinar_info <- sample_info %>%
  dplyr::filter(celltype == "acinar")

gamma_info <- sample_info %>%
  dplyr::filter(celltype == "gamma")

macrophage_info <- sample_info %>%
  dplyr::filter(celltype == "macrophage")

acinar_counts <- counts[acinar_info$id, ]
gamma_counts <- counts[gamma_info$id, ]
macrophage_counts <- counts[macrophage_info$id, ]

acinar_means <- Matrix::colMeans(acinar_counts)
gamma_means <- Matrix::colMeans(gamma_counts)
macrophage_means <- Matrix::colMeans(macrophage_counts)

means_df <- data.frame(
  gene = names(acinar_means),
  acinar_expr = log1p(unname(acinar_means)),
  gamma_expr = log1p(unname(gamma_means)),
  macrophage_expr = log1p(unname(macrophage_means))
)

normalize_bars <- function(LL, FF) {
  
  max_col <- apply(LL, 2, max)
  sweep(FF, 2, max_col, FUN = "*")
  
}

K <- 9
fit <- res_list$pancreas$`1`
colnames(fit$LL) <- paste0(
  "k", c(1, 3, 6, 5, 2, 7, 8, 9, 4)
)
colnames(fit$FF) <- paste0(
  "k", c(1, 3, 6, 5, 2, 7, 8, 9, 4)
)
fit$LL <- fit$LL[,paste0("k", 1:K)]
fit$FF <- fit$FF[,paste0("k", 1:K)]

means_df <- means_df %>%
  dplyr::mutate(
    log1p_k6 = normalize_bars(fit$LL, fit$FF)[, "k6"]
  )

genes_to_label <- names(sort(fit$FF[, "k6"], decreasing = TRUE))[1:10]

g_ag <-  ggplot(data = means_df, aes(x = acinar_expr, y = gamma_expr)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Acinar Cells)") +
  ylab("log(1 + Mean Expression in Gamma Cells)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 4,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.25
  ) + ggtitle("Acinar Vs. Gamma Expression") + cowplot::theme_cowplot()

g_am <-  ggplot(data = means_df, aes(x = acinar_expr, y = macrophage_expr)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Acinar Cells)") +
  ylab("log(1 + Mean Expression in Macrophage Cells)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 4,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.25
  ) + ggtitle("Acinar Vs. Macrophage Expression") + cowplot::theme_cowplot()

g_gm <- ggplot(data = means_df, aes(x = gamma_expr, y = macrophage_expr)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Gamma Cells)") +
  ylab("log(1 + Mean Expression in Macrophage Cells)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 4,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.25
  ) + ggtitle("Gamma Vs. Macrophage Expression") + cowplot::theme_cowplot()

g_fa <- ggplot(data = means_df, aes(x = acinar_expr, y = log1p_k6)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Acinar Cells)") +
  ylab("Factor 6 Score") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 4,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.25
  ) + ggtitle("Acinar Expression Vs. k6") + cowplot::theme_cowplot()

g_fg <- ggplot(data = means_df, aes(x = gamma_expr, y = log1p_k6)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Gamma Cells)") +
  ylab("Factor 6 Score") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 4,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.25
  ) + ggtitle("Gamma Expression Vs. k6") + cowplot::theme_cowplot()

g_fm <- ggplot(data = means_df, aes(x = macrophage_expr, y = log1p_k6)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Macrophage Cells)") +
  ylab("Factor 6 Score") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 4,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.25
  ) + ggtitle("Macrophage Expression Vs. k6") + cowplot::theme_cowplot()


library(ggpubr)

g <- ggarrange(
  g_ag, g_am, g_gm,
  g_fa, g_fm, g_fg,
  nrow = 2, ncol = 3
)

ggsave(
  "~/Downloads/factor_6.pdf",
  g,
  device = "pdf",
  width = 13,
  height = 9
)
