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
colnames(F_tm) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
d_tm <- apply(L_tm,2,max)
F_tm <- scale_cols(F_tm,d_tm)


tm_df <- data.frame(
  k9 = F_tm[,"k9"],
  k12 = F_tm[,"k12"],
  gene = rownames(F_tm)
)

genes_to_label <- c("Vim")

ggplot(data = tm_df, aes(x = log1p(k9), y = log1p(k12))) +
  geom_point(size = 1) +
  geom_text_repel(
    data = subset(tm_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 4,
    max.overlaps = Inf,
    segment.color = "black",     # draw connecting segment
    segment.size = 0.4,
    min.segment.length = 0
  ) 


k9_distinctive <- tm_df %>%
  dplyr::filter(
    log1p(k9) > 3.5 & log1p(k12) < 0.5
  )

genes_to_label <- k9_distinctive$gene

ggplot(data = tm_df, aes(x = log1p(k9), y = log1p(k12))) +
  geom_point(size = 1) +
  geom_text_repel(
    data = subset(tm_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    segment.color = "black",     # draw connecting segment
    segment.size = 0.4,
    min.segment.length = 0
  ) 


k12_distinctive <- tm_df %>%
  dplyr::filter(
    log1p(k12) > 3.5 & log1p(k9) < 0.5
  )

genes_to_label <- k12_distinctive$gene

ggplot(data = tm_df, aes(x = log1p(k9), y = log1p(k12))) +
  geom_point(size = 1) +
  geom_text_repel(
    data = subset(tm_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    segment.color = "black",     # draw connecting segment
    segment.size = 0.4,
    min.segment.length = 0
  ) 


interesting_genes <- c("Fabp4", "Plvap", "Esm1", "Ccl2")

tm_fitted_interesting <- tm_k13$L %*% t(tm_k13$F[interesting_genes,])

S_log1p <- Matrix::Diagonal(x = log1p_k13$s, names = rownames(counts))
log1p_fitted_interesting <- as.matrix(expm1(
  S_log1p %*% log1p_k13$LL %*% t(log1p_k13$FF[interesting_genes,])
))

counts_interesting <- matrix(
  data = c(
    counts[,"Fabp4"],
    counts[,"Plvap"],
    counts[,"Esm1"],
    counts[,"Ccl2"]
  ),
  nrow = nrow(counts),
  ncol = 4
)

plot(log1p(counts_interesting[,1]), log1p(tm_fitted_interesting[,"Fabp4"]))
plot(log1p(counts_interesting[,1]), log1p(log1p_fitted_interesting[,"Fabp4"]))


plot(log1p(counts_interesting[,2]), log1p(tm_fitted_interesting[,"Plvap"]))
plot(log1p(counts_interesting[,2]), log1p(log1p_fitted_interesting[,"Plvap"]))

plot(
  log1p(counts_interesting[,3]), 
  log1p(tm_fitted_interesting[,"Esm1"]),
  xlab = "log1p(Esm1 Raw Count)",
  ylab = "log1p(Esm1 Topic Model Fitted)"
  )
plot(
  log1p(counts_interesting[,3]),
  log1p(log1p_fitted_interesting[,"Esm1"]),
  xlab = "log1p(Esm1 Raw Count)",
  ylab = "log1p(Esm1 log1p Model Fitted)"
  )

plot(
  log1p(counts_interesting[,4]), 
  log1p(tm_fitted_interesting[,"Ccl2"]),
  xlab = "log1p(Ccl2 Raw Count)",
  ylab = "log1p(Ccl2 Topic Model Fitted)"
  )
plot(
  log1p(counts_interesting[,4]), 
  log1p(log1p_fitted_interesting[,"Ccl2"]),
  xlab = "log1p(Ccl2 Raw Count)",
  ylab = "log1p(Ccl2 log1p Model Fitted)"
  )
