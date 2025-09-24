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
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none") +
  ggtitle("Celltype Associated Factors from log1p Model With c = 1")

sp2 <- structure_plot(
  L[i,paste0("k", other_topics)],
  grouping = clusters[i],gap = 15,n = Inf, 
  colors = other_colors)
p2 <- sp2 +
  labs(y = "Membership") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none") +
  ggtitle("Treatment Associated Factors from log1p Model With c = 1")

sp3 <- structure_plot(L[i,paste0("k", other_topics)],grouping = conditions[i],gap = 30,n = Inf,
                      colors = other_colors)
p3 <- sp3 +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none") +
  ggtitle("Treatment Associated Factors from log1p Model With c = 1")


tm_k13 <- res_list$pancreas$`Inf`

L <- poisson2multinom(tm_k13)$L
F <- poisson2multinom(tm_k13)$F
colnames(L) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
L <- L[,paste0("k", 1:13)]

colnames(F) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
F <- F[,paste0("k", 1:13)]


#sp34_loadings_order_call <- structure_plot(
#  L[i,paste0("k", celltype_topics)],grouping = clusters[i])

p4 <- structure_plot(
  L[i,paste0("k", celltype_topics_tm)],grouping = clusters[i],gap = 15
  ) +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor"))	+
  guides(colour = "none") +
  ggtitle("Celltype Associated Factors from Topic Model")				 
p5 <- structure_plot(L[i,paste0("k", other_topics)],grouping = clusters[i],gap = 15,
                     colors = other_colors) +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none") +
  ggtitle("Treatment Associated Factors from Topic Model")
p6 <- structure_plot(L[i,paste0("k", other_topics)],grouping = conditions[i],gap = 30,
                     colors = other_colors) +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none") +
  ggtitle("Treatment Associated Factors from Topic Model")

library(ggpubr)

my_legend <- ggpubr::get_legend(p4)

g1 <- ggarrange(
  p1, p4, common.legend = TRUE, legend = "right",
  legend.grob = my_legend,
  nrow = 2, ncol = 1,
  labels = c("A", "B")
  )

g2 <- ggarrange(
  p2, p5, p3, p6,
  nrow = 4, ncol = 1,
  common.legend = TRUE, legend = "right",
  labels = c("C", "D", "E", "F")
)

g <- ggarrange(g1, g2, heights = c(0.5, 1), nrow = 2, ncol = 1)


ggsave(
  "../images/lsa_structure.png",
  g,
  device = "png",
  width = 11,
  height = 15
)



