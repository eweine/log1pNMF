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

i <- c(sample(which(clusters == "Beta"),900),
       which(clusters != "Beta"))
scale_cols <- function (A, b)
  t(t(A) * b)
celltype_topics <- c(1, 2, 3, 5, 6, 8, 9, 11, 12, 13)
celltype_topics_tm <- c(celltype_topics, 10)
other_topics <- c(4, 7, 10)
L <- log1p_k13$LL
d <- apply(L,2,max)
L <- scale_cols(L,1/d)

colnames(L) <- paste0("k", 1:13)
# swap these two topics for better visibility
l11 <- L[,"k11"]
l1 <- L[,"k1"]
L[,"k11"] <- l1
L[,"k1"] <- l11

tm_k13 <- res_list$pancreas$`Inf`

L2 <- poisson2multinom(tm_k13)$L
colnames(L2) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
L2 <- L2[,paste0("k", 1:13)]

log1p_plot_list <- list()
tm_plot_list <- list()

for (condt in unique(barcodes$condition)) {
  
  # exclude acinar cells because there are so few
  cond_bc <- barcodes %>%
    dplyr::filter(condition == condt & celltype != "Acinar") 
  
  df_nonbeta <- cond_bc %>% filter(celltype != "Beta")
  df_beta    <- cond_bc %>% filter(celltype == "Beta") %>% sample_frac(0.2)
  
  cond_bc <- bind_rows(df_nonbeta, df_beta)
  
  clusters   <- factor(cond_bc$celltype,
                       c("Ductal","Endothelial/Mesnchymal","Macrophage",
                         "Alpha","Beta","Delta","Gamma"))
  
  L_cond <- L[cond_bc$cell_bc, ]
  L2_cond <- L2[cond_bc$cell_bc, ]
  
  log1p_plot_list[[condt]] <- structure_plot(
    L_cond, grouping = clusters, gap = 15, n = Inf
  )
  
  tm_plot_list[[condt]] <- structure_plot(
    L2_cond, grouping = clusters, gap = 15, n = Inf
  )
  
}

library(ggpubr)

g_log1p <- ggarrange(
  plotlist = log1p_plot_list, 
  nrow = 4, 
  ncol = 1, 
  common.legend = TRUE,
  legend = "right",
  labels = names(log1p_plot_list)
  )

g_tm <- ggarrange(
  plotlist = tm_plot_list, 
  nrow = 4, 
  ncol = 1, 
  common.legend = TRUE,
  legend = "right",
  labels = names(log1p_plot_list)
)

ggsave(
  "~/Downloads/lsa_log1p_by_group.png",
  g_log1p,
  device = "png",
  width = 11,
  height = 12
)

ggsave(
  "~/Downloads/lsa_log1p_by_group.png",
  g_tm,
  device = "png",
  width = 11,
  height = 12
)

