library(dplyr)
library(fastTopics)
library(log1pNMF)
library(ggplot2)
library(cowplot)
library(Matrix)
library(readr)
library(stringr)

set.seed(1)

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")
load("../data/experiment_results.Rdata")

s <- Matrix::rowSums(counts)
s <- s / mean(s)

barcodes   <- as.data.frame(barcodes)
barcodes$celltype <- if_else(
  barcodes$celltype == "Endothelial/Mesnchymal",
  "Endothelial/Mesenchymal",
  barcodes$celltype
)

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

i <- c(sample(which(clusters == "Beta"),900),
       which(clusters != "Beta"))
scale_cols <- function (A, b)
  t(t(A) * b)
celltype_topics <- c(1, 2, 3, 5, 6, 8, 9, 11, 12, 13)
celltype_topics_tm <- c(celltype_topics, 10)
other_topics <- c(4, 7, 10)
L <- log1p_k13$LL
FF <- log1p_k13$FF
d <- apply(L,2,max)
L <- scale_cols(L,1/d)

colnames(L) <- paste0("k", 1:13)
colnames(FF) <- paste0("k", 1:13)
# swap these two topics for better visibility
l11 <- L[,"k11"]
l1 <- L[,"k1"]
L[,"k11"] <- l1
L[,"k1"] <- l11

f11 <- FF[,"k11"]
f1 <- FF[,"k1"]
FF[,"k11"] <- f1
FF[,"k1"] <- f11

other_colors <- c("#df8461", "#6f340d", "#00538A")

sp1 <- structure_plot(
  L[i,paste0("k", celltype_topics)],grouping = clusters[i],gap = 15,n = Inf
  )
p1 <- sp1 +
  labs(y = "Membership") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Celltype Associated Factors") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(0.7, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.title.position = "top"
  ) 

sp2 <- structure_plot(
  L[i,paste0("k", other_topics)],
  grouping = clusters[i],gap = 15,n = Inf, 
  colors = other_colors)
p2 <- sp2 +
  labs(y = "Membership") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Treatment Associated Factors") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(1, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.title.position = "top"
  )

sp3 <- structure_plot(L[i,paste0("k", other_topics)],grouping = conditions[i],gap = 30,n = Inf,
                      colors = other_colors)
p3 <- sp3 +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Treatment Associated Factors") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(1, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.title.position = "top"
  )


tm_k13 <- res_list$pancreas$`Inf`

L <- tm_k13$L
FF_tm <- poisson2multinom(tm_k13)$F
colnames(L) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
L <- L[,paste0("k", 1:13)]

colnames(FF_tm) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
FF_tm <- FF_tm[,paste0("k", 1:13)]


#sp34_loadings_order_call <- structure_plot(
#  L[i,paste0("k", celltype_topics)],grouping = clusters[i])

p4 <- structure_plot(
  log1pNMF:::normalize_bars(diag(1/s[i]) %*% L[i,paste0("k", celltype_topics_tm)]),
  grouping = clusters[i],gap = 15
  ) +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor", ncol=1))	+
  guides(colour = "none") +
  ggtitle("c = \u221E Cell Scores - CellType Associated Factors")	+
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(0.7, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.title.position = "top"
  )			 
p5 <- structure_plot(
  log1pNMF:::normalize_bars(diag(1/s[i]) %*% L[i,paste0("k", other_topics)]),
  grouping = clusters[i],gap = 15,
  colors = other_colors
  ) +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none") +
  ggtitle("c = \u221E Cell Scores - Treatment Associated Factors") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(1, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.title.position = "top"
  )
p6 <- structure_plot(
  log1pNMF:::normalize_bars(diag(1/s[i]) %*% L[i,paste0("k", other_topics)]),
  grouping = conditions[i],gap = 30,
  colors = other_colors
  ) +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none") +
  ggtitle("c = \u221E Cell Scores - Treatment Associated Factors") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(1, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.title.position = "top"
  )

library(ggpubr)

my_legend <- ggpubr::get_legend(p4)

g1 <- ggarrange(
  p1, p4, common.legend = TRUE, legend = "right",
  legend.grob = my_legend,
  nrow = 2, ncol = 1,
  labels = c("A", "B")
  )

g2 <- ggarrange(
  p3, p6, p2, p5,
  nrow = 4, ncol = 1,
  common.legend = TRUE, legend = "right",
  labels = c("A", "B", "C", "D")
)

#g <- ggarrange(g1, g2, heights = c(0.5, 1), nrow = 2, ncol = 1)


ggsave(
  "../images/lsa_structure_celltype.png",
  g1,
  device = "png",
  width = 11,
  height = 7
)

ggsave(
  "../images/lsa_structure_treatment.png",
  g2,
  device = "png",
  width = 11,
  height = 14
)

# here, want to get top genes

get_top_genes <- function(f, n_top = 10) {
  
  genes <- names(sort(f, decreasing = TRUE))[1:n_top]
  pasted_genes <- paste(genes, collapse = ", ")
  return(pasted_genes)
  
}

top_genes_log1p <- apply(FF, 2, get_top_genes)
top_genes_log1p <- as.data.frame(top_genes_log1p)
colnames(top_genes_log1p) <- c("top_genes_log1p")
top_genes_log1p$factor <- rownames(top_genes_log1p)
rownames(top_genes_log1p) <- NULL
top_genes_log1p <- top_genes_log1p %>% dplyr::select(c("factor", "top_genes_log1p"))

top_genes_tm <- apply(FF_tm, 2, get_top_genes)
top_genes_tm <- as.data.frame(top_genes_tm)
colnames(top_genes_tm) <- c("top_genes_tm")
top_genes_tm$factor <- rownames(top_genes_tm)
rownames(top_genes_tm) <- NULL
top_genes_tm <- top_genes_tm %>% dplyr::select(c("factor", "top_genes_tm"))

colnames(top_genes_tm) <- c("Factor", "Top Genes - Topic Model")
colnames(top_genes_log1p) <- c("Factor", "Top Genes - log1p Model c = 1")

top_genes <- top_genes_log1p %>%
  dplyr::inner_join(top_genes_tm)

top_genes <- top_genes %>%
  mutate(
    `Top Genes - Topic Model` = str_wrap(`Top Genes - Topic Model`, width = 40),
    `Top Genes - log1p Model c = 1` = str_wrap(`Top Genes - log1p Model c = 1`, width = 40)
  )

top_genes_celltype <- top_genes %>%
  dplyr::filter(
    !(Factor %in% c("k4", "k7", "k10"))
  )

top_genes_treatment <- top_genes %>%
  dplyr::filter(
    Factor %in% c("k4", "k7", "k10")
  )

library(readr)
library(dplyr)
library(stringr)
library(knitr)
library(kableExtra)

df_checker_celltype <- top_genes_celltype %>%
  mutate(
    Factor = Factor,
    # For odd rows, shade only TopWords1
    `Top Genes - Topic Model` = if_else(row_number() %% 2 == 1,
                                        cell_spec(`Top Genes - Topic Model`, "latex", background = "gray!10"),
                                        `Top Genes - Topic Model`
    ),
    # For even rows, shade only TopWords2
    `Top Genes - log1p Model c = 1` = if_else(row_number() %% 2 == 0,
                                                  cell_spec(`Top Genes - log1p Model c = 1`, "latex", background = "gray!10"),
                                                  `Top Genes - log1p Model c = 1`
    )
  )

df_checker_treatment <- top_genes_treatment %>%
  mutate(
    Factor = Factor,
    # For odd rows, shade only TopWords1
    `Top Genes - Topic Model` = if_else(row_number() %% 2 == 1,
                                        cell_spec(`Top Genes - Topic Model`, "latex", background = "gray!10"),
                                        `Top Genes - Topic Model`
    ),
    # For even rows, shade only TopWords2
    `Top Genes - log1p Model c = 1` = if_else(row_number() %% 2 == 0,
                                              cell_spec(`Top Genes - log1p Model c = 1`, "latex", background = "gray!10"),
                                              `Top Genes - log1p Model c = 1`
    )
  )

# 3. Make the LaTeX table
# - 'booktabs = TRUE' for better looking horizontal rules
# - 'escape = FALSE' avoids escaping characters like underscores
# - column_spec can set the width for each column to force wrapping
kable(df_checker_celltype, format = "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(
    position = "center"
  ) %>%
  column_spec(2, width = "5cm") %>%
  column_spec(3, width = "5cm")

kable(df_checker_treatment, format = "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(
    position = "center"
  ) %>%
  column_spec(2, width = "5cm") %>%
  column_spec(3, width = "5cm")

