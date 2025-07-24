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
topic_colors <- c("tomato","darkblue","dodgerblue")
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
  ggtitle("log1p Model Loadings (c = 1)") +
  ylab("Membership") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none")


L0 <- poisson2multinom(tm3_r1_init)$L
L <- L0
L0[,2] <- L[,3]
L0[,3] <- L[,2]
g2 <- structure_plot(L0,grouping = samples$label,topics = rev(1:3),
                     loadings_order = 1:n,colors = topic_colors) +
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 9),
    plot.title = element_text(size = 11)) + 
  ggtitle("Topic Model Loadings") +
  ylab("Membership") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none")

FF_log1p <- log1p_fit3$FF
col_maxima <- apply(log1p_fit3$LL, 2, max)
FF_log1p <- sweep(FF_log1p, 2, col_maxima, FUN = "*")

log1p_df <- data.frame(
  k2 = FF_log1p[,"k2"],
  k3 = FF_log1p[,"k3"]
)

g3 <- ggplot(data = log1p_df, aes(x = k2, y = k3)) +
  geom_point(size = 0.5, alpha = 0.25) +
  xlab("log1p Model k2") +
  ylab("log1p Model k3") +
  theme_cowplot(font_size = 10) +
  ggtitle("log1p Model Factors") +
  theme(plot.title = element_text(hjust = 0.5))


F0 <- tm3_r1_init$F
FF <- F0
F0[,2] <- FF[,3]
F0[,3] <- FF[,2]

L0 <- tm3_r1_init$L
LL <- L0
L0[,2] <- LL[,3]
L0[,3] <- LL[,2]

col_maxima <- apply(L0, 2, max)
F0 <- sweep(F0, 2, col_maxima, FUN = "*")

tm_df <- data.frame(
  k2 = log1p(F0[,2]),
  k3 = log1p(F0[,3])
)

g4 <- ggplot(data = tm_df, aes(x = k2, y = k3)) +
  geom_point(size = 0.5, alpha = 0.25) +
  xlab("log(1 + Topic Model k2)") +
  ylab("log(1 + Topic Model k3)") +
  theme_cowplot(font_size = 10) +
  ggtitle("Topic Model Factors") +
  theme(plot.title = element_text(hjust = 0.5))

g_a <- ggarrange(
  NULL, g0, NULL,
  nrow = 1,
  ncol = 3,
  widths = c(0.5, 1, 0.5),
  labels = c("", "A", "")
)

g_bc <- ggarrange(
  g1, g2,
  nrow = 1,
  ncol = 2,
  common.legend = TRUE,
  legend = "right",
  labels = c("B", "C")
)

g_de <- ggarrange(
  g3, g4,
  nrow = 1,
  ncol = 2,
  labels = c("D", "E")
)

g_final <- ggarrange(
  g_a, g_bc, g_de,
  nrow = 3,
  ncol = 1,
  heights = c(0.8, 0.8, 1)
)

ggsave(
  plot = g_final,
  device = "png",
  filename = "../images/mcf7.png",
  width = 7,
  height = 8
)

