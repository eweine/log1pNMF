library(dplyr)
library(fastTopics)
library(log1pNMF)
library(ggplot2)
library(cowplot)
library(Matrix)
library(readr)

set.seed(1)

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")
load("../data/experiment_results.Rdata")

s <- Matrix::rowSums(counts)
s <- s / mean(s)

barcodes   <- as.data.frame(barcodes)
clusters   <- factor(barcodes$celltype,
                     c("Acinar","Ductal","Endothelial/Mesnchymal","Macrophage",
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

i <- which(clusters == "Endothelial/Mesnchymal")

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
L_norm <-  log1pNMF:::normalize_bars(diag(1/s) %*% L)

#sp34_loadings_order_call <- structure_plot(
#  L[i,paste0("k", celltype_topics)],grouping = clusters[i])

p4 <- structure_plot(
  L_norm[i,paste0("k", celltype_topics_tm)]
)




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
  L[i,paste0("k", celltype_topics)],gap = 15, loadings_order = p4$loadings_order
)

library(ggpubr)

em_plot <- ggarrange(
  sp1$plot,
  p4$plot,
  nrow = 2,
  ncol = 1,
  common.legend = TRUE,
  legend = "right"
)

