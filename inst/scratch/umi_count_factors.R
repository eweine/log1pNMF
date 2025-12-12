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

high11 <- which(L[,"k11"] > 0.25)
high2 <- which(L[,"k2"] > 0.25)

high_both <- intersect(high11, high2)

s <- Matrix::rowSums(counts)

df <- data.frame(
  s = c(
    s[high_both], s[-high_both]
  ),
  class = c(
    rep("k2,k11 high", length(high_both)),
    rep("other", nrow(counts) - length(high_both))
  )
)

ggplot(data =df) +
  geom_density(
    aes(x = s, color = class, fill = class), alpha = 0.5
  ) +
  xlab("Total UMI count")


#counts_high_both <- counts[high_both, ]
#counts_not_high <- counts[-high_both, ]

# make histogram to separate these two out and plot by color

plot(density(Matrix::rowSums(counts_high_both)))


plot(density(Matrix::rowSums(counts)))


