library(rsvd)
library(fastglmpca)
library(fastTopics)
library(flashier)
library(data.table)
library(GEOquery)
library(ggplot2)
library(cowplot)

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

############ Rank 3 Topic Model ############

set.seed(10)
tm3 <- fit_poisson_nmf(counts,k = 3,init.method = "random",
                       numiter = 100,
                       control = list(nc = 7))

topic_colors <- c("tomato","darkblue","dodgerblue")
n <- nrow(counts)
L <- poisson2multinom(tm3)$L
structure_plot(L,grouping = samples$label,topics = 1:3,
               loadings_order = 1:n,colors = topic_colors)$plot +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5))


tm0 <- fastTopics:::fit_pnmf_rank1(counts)

init_LL <- cbind(
  tm0$L,
  matrix(data = 1e-8,
         nrow = nrow(counts),
         ncol = 2)
)

rownames(init_LL) <- rownames(counts)

init_FF <- cbind(
  tm0$F,
  matrix(data = 1e-8,
         nrow = ncol(counts),
         ncol = 2)
)

rownames(init_FF) <- colnames(counts)

tm3_fit0 <- init_poisson_nmf(X = counts, F = init_FF, L = init_LL)
tm3_r1_init <- fit_poisson_nmf(X = counts, fit0 = tm3_fit0, control = list(nc = 7))


topic_colors <- c("tomato","darkblue","dodgerblue")
n <- nrow(counts)
L <- poisson2multinom(tm3_r1_init)$L
structure_plot(L,grouping = samples$label,topics = 1:3,
               loadings_order = 1:n,colors = topic_colors)$plot +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5))

############# Rank 3 Topic model with Rank 1 Initialization ###########

tm0 <- fastTopics:::fit_pnmf_rank1(counts)

init_LL <- cbind(
  tm0$L,
  matrix(data = 1e-8,
         nrow = nrow(counts),
         ncol = 2)
)

rownames(init_LL) <- rownames(counts)

init_FF <- cbind(
  tm0$F,
  matrix(data = 1e-8,
         nrow = ncol(counts),
         ncol = 2)
)

rownames(init_FF) <- colnames(counts)

tm3_fit0 <- init_poisson_nmf(X = counts, F = init_FF, L = init_LL)
tm3_r1_init <- fit_poisson_nmf(X = counts, fit0 = tm3_fit0, control = list(nc = 7, numiter = 1000))


topic_colors <- c("tomato","darkblue","dodgerblue")
n <- nrow(counts)
L <- poisson2multinom(tm3_r1_init)$L
structure_plot(L,grouping = samples$label,topics = 1:3,
               loadings_order = 1:n,colors = topic_colors)$plot +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5))


############# Rank 4 Topic model with Rank 1 Initialization ###########

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
tm4_r1_init <- fit_poisson_nmf(X = counts, fit0 = tm4_fit0, control = list(nc = 7))


topic_colors <- c("olivedrab", "tomato","darkblue","dodgerblue")
n <- nrow(counts)
L <- poisson2multinom(tm4_r1_init)$L
structure_plot(L,grouping = samples$label,topics = 1:4,
               loadings_order = 1:n,colors = topic_colors)$plot +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5))

# the topic model is a bit hard to understand with 4 topics...

############# Rank 3 log1p model with Rank 1 Initialization ###########
library(log1pNMF)

set.seed(10)
log1p_fit3 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 3,
  loglik = "exact",
  control = list(maxiter = 250)
)

normalized_structure_plot(
  log1p_fit3,
  grouping  = samples$label
  )


############# Rank 3 log1p model with random Initialization ###########
set.seed(10)
log1p_fit3 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 3,
  loglik = "exact",
  control = list(maxiter = 250),
  init_method = "random"
)

normalized_structure_plot(
  log1p_fit3,
  grouping  = samples$label
)

############# Rank 4 log1p model with rank1 Initialization ###########

set.seed(10)
log1p_fit4 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 4,
  loglik = "exact",
  control = list(maxiter = 250)
)

normalized_structure_plot(
  log1p_fit4,
  grouping  = samples$label
)

########### Additional Plots #########

# first, compare the two relevant factors for log1p K = 4
FF_log1p <- log1p_fit4$FF
col_maxima <- apply(log1p_fit4$LL, 2, max)
FF_log1p <- sweep(FF_log1p, 2, col_maxima, FUN = "*")

log1p_df <- data.frame(
  k2 = FF_log1p[,"k_2"],
  k3 = FF_log1p[,"k_3"]
)

library(ggplot2)
g1 <- ggplot(data = log1p_df, aes(x = k2, y = k3)) +
  geom_point(size = 0.5, alpha = 0.25) +
  xlab("log1p Model k2") +
  ylab("log1p Model k3")

FF_tm <- tm4_r1_init$F
col_maxima <- apply(tm4_r1_init$L, 2, max)
FF_tm <- sweep(FF_tm, 2, col_maxima, FUN = "*")

tm_df <- data.frame(
  k4 = log1p(FF_tm[,4]),
  k3 = log1p(FF_tm[,3])
)

g2 <- ggplot(data = tm_df, aes(x = k4, y = k3)) +
  geom_point(size = 0.5, alpha = 0.25) +
  xlab("log1p(Topic Model k4)") +
  ylab("log1p(Topic Model k3)")

library(ggpubr)
ggarrange(g1, g2, nrow = 1, ncol = 2)


inter_RA_df <- data.frame(
  log1p_k2 = FF_log1p[,"k_3"],
  tm_k4 = log1p(FF_tm[,3])
)

ggplot(data = inter_RA_df, aes(x = log1p_k2, y = tm_k4)) +
  geom_point(size = 0.5, alpha = 0.25) +
  xlab("log1p(Topic Model k4)") +
  ylab("log1p(Topic Model k3)")


tm_df2 <- data.frame(
  ra_tgfb = log1p(FF_tm[,2]),
  tgfb = log1p(FF_tm[,3]),
  ra = log1p(FF_tm[,4]),
  ra_tgfb_avg = log1p(0.5 * FF_tm[,3] + 0.5 * FF_tm[,4])
)

ggplot(data = tm_df2, aes(x = ra_tgfb_avg, y = ra_tgfb)) +
  geom_point(size = 0.5, alpha = 0.25) +
  xlab("log1p(Topic Model TGFB)") +
  ylab("log1p(Topic Model RA+TGFB)")


