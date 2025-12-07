library(rsvd)
library(fastglmpca)
library(fastTopics)
library(flashier)
library(data.table)
library(GEOquery)
library(ggplot2)
library(cowplot)
library(log1pNMF)
library(ggpubr)
library(DESeq2)

load("../data/experiment_results.Rdata")
set.seed(10)
counts <- fread("../data/raw_data/GSE152749_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(counts) <- "data.frame"
rownames(counts) <- counts$GeneID
counts <- counts[,-1]
counts <- as.matrix(counts)
storage.mode(counts) <- "double"
counts <- t(counts)
ids <- rownames(counts)

genes <- fread("../data/raw_data/Human.GRCh38.p13.annot.tsv.gz",sep = "\t",
               header = TRUE,stringsAsFactors = FALSE)
class(genes) <- "data.frame"
genes <- genes[1:10]
genes <- transform(genes,
                   GeneType = factor(GeneType),
                   Status   = factor(Status))

geo <- getGEO(filename = "../data/raw_data/GSE152749_family.soft.gz")
samples <- data.frame(id = names(GSMList(geo)),
                      treatment = sapply(GSMList(geo),
                                         function (x) Meta(x)$title))
samples <- samples[ids,]
rownames(samples) <- NULL
samples <- transform(samples,
                     EtOH = grepl("EtOH",treatment,fixed = TRUE),
                     RA   = grepl("RA",treatment,fixed = TRUE),
                     TGFb = grepl("TGFb",treatment,fixed = TRUE))
samples$label                 <- "EtOH"
samples[samples$RA,"label"]   <- "RA"
samples[samples$TGFb,"label"] <- "TGFb"
samples[with(samples,RA & TGFb),"label"] <- "RA+TGFb"
samples <- transform(samples,
                     label = factor(label,c("EtOH","RA","TGFb","RA+TGFb")))


x <- colSums(counts > 0)
i <- which(x > 3 &
             genes$GeneType == "protein-coding" &
             genes$Status == "active")
genes  <- genes[i,]
counts <- counts[,i]
colnames(counts) <- genes$Symbol

# Run DE-Seq
# ------------------------------------------------------------
# 1. Construct DESeq2 dataset
# ------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = t(counts),   # DESeq2 expects genes × samples
  colData   = samples,
  design    = ~ label
)

# ------------------------------------------------------------
# 2. Run DESeq2 pipeline
# ------------------------------------------------------------
dds <- DESeq(dds)

# ------------------------------------------------------------
# 3. Get the contrasts
# ------------------------------------------------------------
# Contrast syntax: c("label", "level_of_interest", "reference")

# RA vs EtOH
res_RA_vs_EtOH <- as.data.frame(results(dds, contrast = c("label", "RA", "EtOH")))

# TGFb vs EtOH
res_TGFb_vs_EtOH <- as.data.frame(results(dds, contrast = c("label", "TGFb", "EtOH")))




set.seed(1)
fgpca_fit <- res_list$mcf7$glmpca

pdat <- data.frame(samples,fgpca_fit$V)
pdat$Group <- pdat$label
g0 <- ggplot(pdat,aes(x = k_1,y = k_2,color = Group)) +
  geom_point() +
  scale_color_manual(values = c("tomato","darkblue","dodgerblue",
                                "olivedrab")) + 
  theme_cowplot(font_size = 10) +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("GLM-PCA Loadings") +
  theme(
    plot.title = element_text(hjust = 0.5)
    )

log1p_fit3 <- res_list$mcf7$`1`

tm3_r1_init <- res_list$mcf7$`Inf`

n <- nrow(counts)
colnames(log1p_fit3$FF) <- paste0("k", 1:3)
colnames(log1p_fit3$LL) <- paste0("k", 1:3)
topic_colors <- c("black","orange","blue")
g1 <- normalized_structure_plot(
  log1p_fit3,
  grouping = samples$label,
  loadings_order = 1:n,
  colors = topic_colors,
  topics = rev(1:3)
) + 
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 9),
    plot.title = element_text(size = 11)) + 
  ggtitle("c = 1 Sample Scores") +
  ylab("Membership") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none")


L0 <- tm3_r1_init$L
L <- L0
L0[,2] <- L[,3]
L0[,3] <- L[,2]
g2 <- structure_plot(
  log1pNMF:::normalize_bars(diag(1 / log1p_fit3$s) %*% L0),
  grouping = samples$label,topics = rev(1:3),
  loadings_order = 1:n,colors = topic_colors
  ) +
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 9),
    plot.title = element_text(size = 11)) + 
  ggtitle("c = \u221E Sample Scores") +
  ylab("Membership") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none")

FF_log1p <- log1p_fit3$FF
col_maxima <- apply(log1p_fit3$LL, 2, max)
FF_log1p <- sweep(FF_log1p, 2, col_maxima, FUN = "*")

log1p_df <- data.frame(
  k1 = FF_log1p[,"k1"],
  k2 = FF_log1p[,"k2"],
  k3 = FF_log1p[,"k3"],
  lfc_ra = res_RA_vs_EtOH$log2FoldChange,
  padj_ra = res_RA_vs_EtOH$padj,
  lfc_tgfb = res_TGFb_vs_EtOH$log2FoldChange,
  padj_tgfb = res_TGFb_vs_EtOH$padj
)

log1p_df$col <- "lightgrey"   # default

# RA upregulated
log1p_df$col[log1p_df$padj_ra < 0.01 & log1p_df$lfc_ra > 1] <- "orange"

# TGFb upregulated
log1p_df$col[log1p_df$padj_tgfb < 0.01 & log1p_df$lfc_tgfb > 1] <- "blue"

# Both upregulated
log1p_df$col[
  log1p_df$padj_tgfb < 0.01 & log1p_df$lfc_tgfb > 1 &
    log1p_df$padj_ra < 0.01 & log1p_df$lfc_ra > 1
  ] <- "forestgreen"

g3 <- ggplot(data = log1p_df, aes(x = k2, y = k3)) +
  geom_point(aes(color = col), size = 0.5, alpha = 2/3) +
  scale_color_identity() +
  xlab("k2") +
  ylab("k3") +
  theme_cowplot(font_size = 10) +
  ggtitle("c = 1 Gene Scores") +
  theme(plot.title = element_text(hjust = 0.5))


F0 <- tm3_r1_init$F
FF <- F0
F0[,2] <- FF[,3]
F0[,3] <- FF[,2]

L0 <- tm3_r1_init$L
LL <- L0
L0[,2] <- LL[,3]
L0[,3] <- LL[,2]
L0_norm <- diag(1 / log1p_fit3$s) %*% L0

col_maxima <- apply(L0, 2, max)
F0 <- sweep(F0, 2, col_maxima, FUN = "*")

tm_df <- data.frame(
  k1 = F0[,1],
  k2 = F0[,2],
  k3 = F0[,3],
  lfc_ra = res_RA_vs_EtOH$log2FoldChange,
  padj_ra = res_RA_vs_EtOH$padj,
  lfc_tgfb = res_TGFb_vs_EtOH$log2FoldChange,
  padj_tgfb = res_TGFb_vs_EtOH$padj
)

rownames(F0) <- genes$Symbol
rownames(FF_log1p) <- genes$Symbol

tm_df$col <- "lightgrey"   # default

# RA upregulated
tm_df$col[tm_df$padj_ra < 0.01 & tm_df$lfc_ra > 1] <- "orange"

# TGFb upregulated
tm_df$col[tm_df$padj_tgfb < 0.01 & tm_df$lfc_tgfb > 1] <- "blue"

# Both upregulated
tm_df$col[
  tm_df$padj_tgfb < 0.01 & tm_df$lfc_tgfb > 1 &
    tm_df$padj_ra < 0.01 & tm_df$lfc_ra > 1
] <- "forestgreen"

g4 <- ggplot(data = tm_df, aes(x = k2, y = k3)) +
  geom_point(aes(color = col), size = 0.5, alpha = 2/3) +
  scale_color_identity() +
  xlab("k2") +
  ylab("k3") +
  theme_cowplot(font_size = 10) +
  ggtitle("c = \u221E Gene Scores") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 100, 1000, 10000)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 100, 1000, 10000)) +
  theme(plot.title = element_text(hjust = 0.5)) 

g_bc <- ggarrange(
  g2, g1,
  nrow = 1,
  ncol = 2,
  common.legend = TRUE,
  legend = "right",
  labels = c("A", "B")
)

g_de <- ggarrange(
  g4, g3,
  nrow = 1,
  ncol = 2,
  labels = c("C", "D")
)

g_final <- ggarrange(
  g_bc, g_de,
  nrow = 2,
  ncol = 1,
  heights = c(0.8, 1)
)

ggsave(
  plot = g_final,
  device = "png",
  filename = "../images/mcf7.png",
  width = 7,
  height = 6
)

# additional exploration regarding factors here

log1p_df$gene <- genes$Symbol
tm_df$gene <- genes$Symbol

# top 4 genes in log1p k2 are known to interact / metabolize with retinoic acid
#tm_df$diff <- tm_df$k2 - tm_df$k3
#tm_df$rel_diff <- log1p(tm_df$k2) - log1p(tm_df$k3)

top_genes <- data.frame(
  factor = c("k1", "k2", "k3"),
  top_genes = c(
    paste(
      log1p_df %>% dplyr::arrange(desc(k1)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(gene),
      collapse = ", "
    ),
    paste(
      log1p_df %>% dplyr::arrange(desc(k2)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(gene),
      collapse = ", "
    ),
    paste(
      log1p_df %>% dplyr::arrange(desc(k3)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(gene),
      collapse = ", "
    )
  )
)

top_genes_tm <- data.frame(
  factor = c("k1", "k2", "k3"),
  top_genes = c(
    paste(
      tm_df %>% dplyr::arrange(desc(k1)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(gene),
      collapse = ", "
    ),
    paste(
      tm_df %>% dplyr::arrange(desc(k2)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(gene),
      collapse = ", "
    ),
    paste(
      tm_df %>% dplyr::arrange(desc(k3)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(gene),
      collapse = ", "
    )
  )
)

colnames(top_genes_tm) <- c("Factor", "Top Genes - Topic Model")
colnames(top_genes) <- c("Factor", "Top Genes - log1p Model c = 1")

top_genes <- top_genes_tm %>%
  dplyr::inner_join(top_genes)



library(kableExtra)

library(readr)
library(dplyr)
library(stringr)
library(knitr)
library(kableExtra)

top_genes <- top_genes %>%
  mutate(
    `Top Genes - Topic Model` = str_wrap(`Top Genes - Topic Model`, width = 40),
    `Top Genes - log1p Model c = 1` = str_wrap(`Top Genes - log1p Model c = 1`, width = 40)
  )

df_checker <- top_genes %>%
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
kable(df_checker, format = "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(
    position = "center"
  ) %>%
  column_spec(2, width = "5cm") %>%
  column_spec(3, width = "5cm")

