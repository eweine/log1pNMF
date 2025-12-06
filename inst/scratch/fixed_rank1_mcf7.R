library(rsvd)
library(fastglmpca)
library(fastTopics)
library(flashier)
library(data.table)
library(GEOquery)
library(ggplot2)
library(cowplot)
library(log1pNMF)

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

n <- nrow(counts)
p <- ncol(counts)

set.seed(1)
LL <- matrix(
  data = runif(n = n * 2, min = 0, max = 0.01),
  nrow = n,
  ncol = 2
)
LL <- cbind(rep(1, n), LL)

FF <- matrix(
  data = runif(n = p * 2, min = 0, max = 0.01),
  nrow = p,
  ncol = 2
)
FF <- cbind(log1p(colMeans(counts)), FF)


log1p_mod_fixed_factor <- fit_poisson_log1p_nmf(
  Y = counts, loglik = "exact",
  init_LL = LL, init_FF = FF,
  update_idx_LL = 2:3,
  update_idx_FF = 1:3,
  control = list(maxiter = 2500, verbose = TRUE),
  s = rep(1, nrow(counts))
)

rownames(log1p_mod_fixed_factor$FF) <- genes$Symbol

log1p_mod_fixed_factor$LL <- log1p_mod_fixed_factor$LL %*% diag(1/colMeans(log1p_mod_fixed_factor$LL))
log1p_mod_fixed_factor$FF <- log1p_mod_fixed_factor$FF %*% diag(colMeans(log1p_mod_fixed_factor$LL))

plot(log1p_mod_fixed_factor$FF[,1], log1p_mod_fixed_factor$FF[,2])

df <- data.frame(
  factor = c(
    rep("k1", nrow(log1p_mod_fixed_factor$FF)),
    rep("k2", nrow(log1p_mod_fixed_factor$FF))
  ),
  value = c(
    log1p_mod_fixed_factor$FF[,1],
    log1p_mod_fixed_factor$FF[,2]
  )
)

ggplot(data = df) +
  geom_density(aes(x = value, color = factor))


df2 <- data.frame(
  k1 = log1p_mod_fixed_factor$FF[,1],
  k2 = log1p_mod_fixed_factor$FF[,2]
)

ggplot(data = df2) +
  geom_point(aes(x = k1, y = k2), size = 0.5, alpha = 0.5)

set.seed(10)
log1p_mod_std <- fit_poisson_log1p_nmf(
  Y = counts, loglik = "exact",
  init_method = "rank1",
  K = 3,
  control = list(maxiter = 1000, verbose = TRUE),
  s = rep(1, nrow(counts))
)

topic_colors <- c("tomato","darkblue","dodgerblue")

normalized_structure_plot(
  log1p_mod_fixed_factor,
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

normalized_structure_plot(
  log1p_mod_std,
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


rownames(log1p_mod_fixed_factor$FF) <- genes$Symbol


