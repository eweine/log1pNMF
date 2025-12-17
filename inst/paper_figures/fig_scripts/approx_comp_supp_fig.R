library(dplyr)
library(fastTopics)
library(log1pNMF)
library(ggplot2)
library(cowplot)
library(Matrix)
library(readr)
library(ggrepel)
library(ggpubr)

hoyer_sparsity <- function(x) {
  
  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))
  
}

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

sparsity_exact_L <- apply(
  L, 2, hoyer_sparsity
)

sparsity_exact_F <- apply(
  FF, 2, hoyer_sparsity
)

cor_exact_F <- cor(FF, method = "spearman")

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
  ggtitle("c = 1 Cell Scores - Celltype Associated Factors (Fitted w/ Exact Log-Likelihood)") +
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
  grouping = conditions[i],gap = 30,n = Inf,
  colors = other_colors,
  topics = paste0(
    "k",
    c(10, 4, 7)
  )
)
p2 <- sp2 +
  labs(y = "Membership",fill = "") +
  guides(fill="none") +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Treatment Associated Factors (Fitted w/ Exact Log-Likelihood)") +
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
  c(11,2,13,3,5,6,7,9,8,10,1,12,4)
)
L <- L[,paste0("k", 1:13)]
d <- apply(L,2,max)
L <- scale_cols(L,1/d)

colnames(FF_log1p_cheb) <- paste0(
  "k", 
  c(11,2,13,3,5,6,7,9,8,10,1,12,4)
)
FF_log1p_cheb <- FF_log1p_cheb[,paste0("k", 1:13)]

sparsity_cheb_L <- apply(
  L, 2, hoyer_sparsity
)

sparsity_cheb_F <- apply(
  FF_log1p_cheb, 2, hoyer_sparsity
)

cor_cheb_F <- cor(FF_log1p_cheb, method = "spearman")

sp3 <- structure_plot(
  L[i, celltype_topics],
  grouping = clusters[i],
  gap = 20,
  topics = rev(paste0(
    "k",
    c(8, 11, 2, 13, 3, 9, 12, 6, 1, 5)
  ))
)

p3 <- sp3 +
  labs(y = "Membership") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Celltype Associated Factors (Fitted w/ Approximate Log-Likelihood)") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(0.7, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.title.position = "top"
  ) 

sp4 <- structure_plot(
  L[i, other_topics],
  grouping = conditions[i],
  gap = 20,
  colors = other_colors,
  topics = paste0(
    "k",
    c(10, 4, 7)
  )
)
p4 <- sp4 +
  labs(y = "Membership",fill = "") +
  guides(fill="none") +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Treatment Associated Factors (Fitted w/ Approximate Log-Likelihood)") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(1, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )

log1p_frob <- res_list$pancreas$`c = 1, frob`

L <- log1p_frob$W
FF_log1p_frob <- t(log1p_frob$H)
colnames(L) <- paste0(
  "k", 
  c(11,6,2,4,5,12,1,3,8,10,13,9,7)
)
L <- L[,paste0("k", 1:13)]
d <- apply(L,2,max)
L <- scale_cols(L,1/d)

colnames(FF_log1p_frob) <- paste0(
  "k", 
  c(11,6,2,4,5,12,1,3,8,10,13,9,7)
)
FF_log1p_frob <- FF_log1p_frob[,paste0("k", 1:13)]

other_topics_frob <- c(
  4, 7, 10, 13
)

sparsity_frob_L <- apply(
  L, 2, hoyer_sparsity
)

sparsity_frob_F <- apply(
  FF_log1p_frob, 2, hoyer_sparsity
)

cor_frob_F <- cor(FF_log1p_frob, method = "spearman")

sp5 <- structure_plot(
  L[i, celltype_topics],
  grouping = clusters[i],
  gap = 20
)
p5 <- sp5 +
  labs(y = "Membership") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Celltype Associated Factors (Fitted w/ Frob NMF on log1p Transformed Counts)") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(0.7, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.title.position = "top"
  ) 


sp6 <- structure_plot(
  L[i, other_topics_frob],
  grouping = conditions[i],
  gap = 20,
  colors = c("#df8461", "#6f340d", "#00538A", "#F6768E"),
  topics = paste0(
    "k",
    c(13, 10, 4, 7)
  )
)

p6 <- sp6 +
  labs(y = "Membership",fill = "") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none") +
  ggtitle("c = 1 Cell Scores - Treatment Associated Factors (Fitted w/ Frob NMF on log1p Transformed Counts)") +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.key.size = unit(1, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )

g1 <- ggarrange(
  p1, p3, p5, 
  nrow = 3, ncol = 1, 
  legend = "right", common.legend = TRUE
  )

g2 <- ggarrange(
  p2, p4, p6, 
  nrow = 3, ncol = 1, 
  legend.grob = get_legend(p6),
  legend = "right", common.legend = TRUE
)

ggsave(
  "../images/supp_approx_lsa_structure_celltype.png",
  g1,
  device = "png",
  width = 11,
  height = 11
)

ggsave(
  "../images/supp_approx_lsa_structure_treatment.png",
  g2,
  device = "png",
  width = 11,
  height = 11
)

# analysis for Frobenius K = 12 fit here:
log1p_frob <- res_list$pancreas$`c = 1, frob, k = 12`

L <- log1p_frob$W
FF_log1p_frob <- t(log1p_frob$H)
colnames(L) <- paste0(
  "k", 
  1:12
)
L <- L[,paste0("k", 1:12)]
d <- apply(L,2,max)
L <- scale_cols(L,1/d)

colnames(FF_log1p_frob) <- paste0(
  "k", 
  1:12
)
FF_log1p_frob <- FF_log1p_frob[,paste0("k", 1:12)]

other_topics <- c(
  4, 9, 8, 11
)

structure_plot(
  L[i, other_topics],
  grouping = conditions[i],
  gap = 20
)

log1p_frob <- res_list$pancreas$`c = 1, frob, k = 11`

L <- log1p_frob$W
FF_log1p_frob <- t(log1p_frob$H)
colnames(L) <- paste0(
  "k", 
  1:11
)
L <- L[,paste0("k", 1:11)]
d <- apply(L,2,max)
L <- scale_cols(L,1/d)

colnames(FF_log1p_frob) <- paste0(
  "k", 
  1:11
)
FF_log1p_frob <- FF_log1p_frob[,paste0("k", 1:11)]

other_topics <- c(
  1, 6, 9, 11
)

structure_plot(
  L[i, other_topics],
  grouping = conditions[i],
  gap = 20
)



