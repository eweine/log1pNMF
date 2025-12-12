# NOTE: all the colors are now aligned properly, I just need to make aesthetic
# updates to the figure

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

barcodes   <- as.data.frame(barcodes)

barcodes <- barcodes %>%
  dplyr::mutate(
    condition = if_else(
      condition == "IL-1B_IFNg",
      "IL-1B + IFNg",
      condition
    )
  )

barcodes$celltype <- if_else(
  barcodes$celltype == "Endothelial/Mesnchymal",
  "Endothelial/Mesenchymal",
  barcodes$celltype
)

clusters   <- factor(barcodes$celltype,
                     c("Acinar","Ductal","Endothelial/Mesenchymal","Macrophage",
                       "Alpha","Beta","Delta","Gamma"))

conditions <- factor(barcodes$condition,
                     c("Untreated","IL-1B","IFNg","IL-1B + IFNg"))

log1p_k13 <- res_list$pancreas$`1`

i <- c(sample(which(clusters == "Beta"),900),
       which(clusters != "Beta"))
scale_cols <- function (A, b)
  t(t(A) * b)
celltype_topics <- c(1, 2, 3, 5, 6, 8, 9, 11, 12, 13)
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
  ggtitle("c = 1 Cell Scores - Celltype Associated Factors") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(0.7, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )

sp2 <- structure_plot(L[i,paste0("k", other_topics)],grouping = conditions[i],gap = 30,n = Inf,
                      colors = other_colors)
p2 <- sp2 +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Treatment Associated Factors") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(1, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )
# now that I have the above, I want to make the same plots for both of
# the approximate fits

log1p_cheb <- res_list$pancreas$`c = 1, cheby`

L <- log1p_cheb$LL
FF_log1p_cheb <- log1p_cheb$FF
colnames(L) <- paste0(
  "k", 
  c(11,2,3,13,5,6,7,9,8,10,1,12,4)
)
L <- L[,paste0("k", 1:13)]
d <- apply(L,2,max)
L <- scale_cols(L,1/d)

colnames(FF_log1p_cheb) <- paste0(
  "k", 
  c(11,2,3,13,5,6,7,9,8,10,1,12,4)
)
FF_log1p_cheb <- FF_log1p_cheb[,paste0("k", 1:13)]


p3 <- structure_plot(
  L[i, celltype_topics],
  grouping = clusters[i],
  gap = 20,
  topics = rev(paste0(
    "k",
    c(8, 11, 2, 13, 3, 9, 12, 6, 1, 5)
  ))
)

p4 <- structure_plot(
  L[i, other_topics],
  grouping = conditions[i],
  gap = 20,
  colors = other_colors
)

log1p_frob <- res_list$pancreas$`c = 1, frob`

L <- log1p_frob$W
FF_log1p_frob <- t(log1p_frob$H)
colnames(L) <- paste0(
  "k", 
  c(11,6,2,4,5,12,1,9,8,7,10,3,13)
)
L <- L[,paste0("k", 1:13)]
d <- apply(L,2,max)
L <- scale_cols(L,1/d)

colnames(FF_log1p_frob) <- paste0(
  "k", 
  c(1,2,3,4,5,6,7,8,9,10,11,12,13)
)
FF_log1p_frob <- FF_log1p_frob[,paste0("k", 1:13)]

other_topics_frob <- c(
  4, 7, 10, 13
)

p5 <- structure_plot(
  L[i, celltype_topics],
  grouping = clusters[i],
  gap = 20
)
p6 <- structure_plot(
  L[i, other_topics_frob],
  grouping = conditions[i],
  gap = 20,
  colors = c(fastTopics:::kelly()[-1][13], "#00538A",  "#df8461", "#6f340d"),
  topics = paste0(
    "k",
    c(13, 10, 4, 7)
  )
)

ggarrange(p5, p3, nrow = 2, ncol = 1, legend = "right")

ggarrange(p2, p6, nrow = 2, ncol = 1, legend = "right")



# these two look great, I just need to add some aesthetic changes
ggarrange(sp1, p3, nrow = 2, ncol = 1, legend = "right")
ggarrange(sp2, p4, nrow = 2, ncol = 1, legend = "right")

