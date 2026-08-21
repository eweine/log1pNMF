library(dplyr)
library(fastTopics)
library(log1pNMF)
library(ggplot2)
library(cowplot)
library(Matrix)
library(readr)
library(ggrepel)
library(ggpubr)

set.seed(1)

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")
load("../data/experiment_results.Rdata")

s <- Matrix::rowSums(counts)
s <- s / mean(s)

barcodes   <- as.data.frame(barcodes)
clusters   <- factor(barcodes$celltype,
                     c("Acinar","Ductal","Endothelial/Mesenchymal","Macrophage",
                       "Alpha","Beta","Delta","Gamma"))

barcodes <- barcodes %>%
  dplyr::mutate(
    condition = if_else(
      condition == "IL-1B_IFNg",
      "IL-1B + IFNg",
      condition
    )
  )

conditions <- factor(barcodes$condition,
                     c("Untreated","IL-1B","IFNg","IL-1B + IFNg"))

log1p_k13 <- res_list$pancreas$`1`

scale_cols <- function (A, b)
  t(t(A) * b)
F_log1p <- log1p_k13$FF
L_log1p <- log1p_k13$LL
d_log1p <- apply(L_log1p,2,max)
F_log1p <- scale_cols(F_log1p,d_log1p)

colnames(F_log1p) <- paste0("k", 1:13)
# swap these two topics for better visibility
f11 <- F_log1p[,"k11"]
f1 <- F_log1p[,"k1"]
F_log1p[,"k11"] <- f1
F_log1p[,"k1"] <- f11

tm_k13 <- res_list$pancreas$`Inf`

F_tm <- tm_k13$F
L_tm <- tm_k13$L
L_tm <- Matrix::Diagonal(x = 1/s) %*% L_tm
## match to the log1p c = 1 reference (same procedure and reference as every
## other pancreas figure; see factor_matching.R). Labels k3/k13 -- the two
## used below -- are identical under the old manual vector and the matcher.
source("factor_matching.R")
tm_labels <- match_display_labels(
  res_list$pancreas[["1"]]$FF, tm_k13$F,
  c(11, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 12, 13))
colnames(F_tm) <- paste0("k", tm_labels)
d_tm <- apply(L_tm,2,max)
F_tm <- scale_cols(F_tm,d_tm)

########## Start new code here ########## 

delta_counts <- counts[which(clusters == "Delta"), ]
gamma_counts <- counts[which(clusters == "Gamma"), ]

gene_means_df <- data.frame(
  gene = colnames(delta_counts),
  delta = unname(Matrix::colMeans(delta_counts)),
  gamma = unname(Matrix::colMeans(gamma_counts)),
  log1p_delta = unname(F_log1p[,"k3"]),
  tm_delta = unname(F_tm[,"k3"]),
  log1p_gamma = unname(F_log1p[,"k13"]),
  tm_gamma = unname(F_tm[,"k13"])
)

genes_to_label <- c("Ppy", "Pyy", "Iapp", "Sst", "Rbp4")

g3 <- ggplot(data = gene_means_df, aes(x = delta, y = gamma)) +
  geom_point(size = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_text_repel(
    data = subset(gene_means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 3.5,
    max.overlaps = Inf,
    segment.color = "black",     # draw connecting segment
    segment.size = 0.4,
    min.segment.length = 0
  ) +
  scale_x_continuous(
    trans = "log1p", 
    breaks = c(0, 10, 500, 2000),
    limits = c(0, 2000)
    ) +
  scale_y_continuous(
    trans = "log1p", 
    breaks = c(0, 10, 500, 2000, 4000),
    limits = c(0, 4500)
    ) +
  xlab("Mean Expression in Delta Cells") +
  ylab("Mean Expression in Gamma Cells") +
  theme_cowplot() +
  ggtitle("Raw Data") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
    )  

g1 <- ggplot(data = gene_means_df, aes(x = log1p_delta, y = log1p_gamma)) +
  geom_point(size = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_text_repel(
    data = subset(gene_means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 3.5,
    max.overlaps = Inf,
    segment.color = "black",     # draw connecting segment
    segment.size = 0.4,
    min.segment.length = 0,
    force_pull = 2,
    nudge_x = 0.5, 
    nudge_y = 0.75
  ) +
  theme_cowplot(font_size = 12) +
  ggtitle("c = 1 Gene Scores") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  ) +
  xlab("k3") +
  ylab("k13") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  )

g2 <- ggplot(data = gene_means_df, aes(x = tm_delta, y = tm_gamma)) +
  geom_point(size = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_text_repel(
    data = subset(gene_means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 3.5,
    max.overlaps = Inf,
    segment.color = "black",     # draw connecting segment
    segment.size = 0.4,
    min.segment.length = 0,
    force_pull = 2,
    nudge_x = 0.5, 
    nudge_y = 0.5
  ) +
  scale_x_continuous(
    trans = "log1p", 
    breaks = c(0, 10, 500, 2000)
  ) +
  scale_y_continuous(
    trans = "log1p", 
    breaks = c(0, 10, 500, 2000, 5000, 10000)
  ) +
  theme_cowplot() +
  ggtitle("c = \u221E Gene Scores") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  ) +
  xlab("k3") +
  ylab("k13")  +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  )

## one-row layout (A, B, C side by side) to halve the figure's page height
g <- ggarrange(
  g1, g2, g3,
  nrow = 1, ncol = 3,
  labels = c("A", "B", "C")
)

ggsave(
  plot = g,
  device = "png",
  filename = "../images/lsa_delta_gamma.png",
  width = 11.5,
  height = 4
)
