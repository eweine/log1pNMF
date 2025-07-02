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

set.seed(10)
counts <- fread("../data/GSE152749_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(counts) <- "data.frame"
rownames(counts) <- counts$GeneID
counts <- counts[,-1]
counts <- as.matrix(counts)
storage.mode(counts) <- "double"
counts <- t(counts)
ids <- rownames(counts)

genes <- fread("../data/Human.GRCh38.p13.annot.tsv.gz",sep = "\t",
               header = TRUE,stringsAsFactors = FALSE)
class(genes) <- "data.frame"
genes <- genes[1:10]
genes <- transform(genes,
                   GeneType = factor(GeneType),
                   Status   = factor(Status))

geo <- getGEO(filename = "../data/GSE152749_family.soft.gz")
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
fgpca_fit <- fit_glmpca_pois(
  Y = t(counts),
  K = 2,
  control = list(maxiter = 250)
)

pdat <- data.frame(samples,fgpca_fit$V)
pdat$Group <- pdat$label
g0 <- ggplot(pdat,aes(x = k_1,y = k_2,color = Group)) +
  geom_point() +
  scale_color_manual(values = c("olivedrab","tomato","darkblue",
                                "dodgerblue")) + 
  theme_cowplot(font_size = 10) +
  xlab("PC1") +
  ylab("PC2") 

set.seed(10)
log1p_fit4 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 4,
  loglik = "exact",
  control = list(maxiter = 250)
)

tm0 <- fastTopics:::fit_pnmf_rank1(counts)

init_LL <- cbind(
  tm0$L,
  matrix(data = 1e-8,
         nrow = nrow(counts),
         ncol = 3)
)

rownames(init_LL) <- rownames(counts)

init_FF <- cbind(
  tm0$F,
  matrix(data = 1e-8,
         nrow = ncol(counts),
         ncol = 3)
)

rownames(init_FF) <- colnames(counts)

tm4_fit0 <- init_poisson_nmf(X = counts, F = init_FF, L = init_LL)
tm4_r1_init <- fit_poisson_nmf(
  X = counts, 
  fit0 = tm4_fit0, 
  numiter = 250,
  control = list(nc = 7)
  )

n <- nrow(counts)
colnames(log1p_fit4$FF) <- paste0("k", 1:4)
colnames(log1p_fit4$LL) <- paste0("k", 1:4)
topic_colors <- c("olivedrab", "tomato","darkblue","dodgerblue")
g1 <- normalized_structure_plot(
  log1p_fit4,
  grouping = samples$label,
  loadings_order = 1:n,
  colors = topic_colors,
  topics = c(2,3,4,1)
)$plot + 
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5)) + 
  ggtitle("log1p Model With c = 1") 


L0 <- poisson2multinom(tm4_r1_init)$L
L <- L0
L0[,2] <- L[,4]
L0[,4] <- L[,2]
g2 <- structure_plot(L0,grouping = samples$label,topics = c(2,3,4,1),
               loadings_order = 1:n,colors = topic_colors)$plot +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5)) + 
  ggtitle("Topic Model") 

FF_log1p <- log1p_fit4$FF
col_maxima <- apply(log1p_fit4$LL, 2, max)
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
  ggtitle("log1p Model Factor Correlation") +
  theme(plot.title = element_text(hjust = 0.5))


F0 <- tm4_r1_init$F
FF <- F0
F0[,2] <- FF[,4]
F0[,4] <- FF[,2]

L0 <- tm4_r1_init$L
LL <- L0
L0[,2] <- LL[,4]
L0[,4] <- LL[,2]

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
  ggtitle("Topic Model Factor Correlation") +
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
  filename = "~/Documents/log1pNMF/inst/paper_figures/pdfs/mcf7.png",
  width = 7,
  height = 8
  )

