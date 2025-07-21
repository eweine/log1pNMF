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
tm3_r1_init <- fit_poisson_nmf(
  X = counts, 
  fit0 = tm3_fit0, 
  numiter = 1000,
  control = list(nc = 7)
)

ll_tm3 <- sum(loglik_poisson_nmf(X = counts, fit = tm3_r1_init))
ll_per_pt_tm3 <- loglik_poisson_nmf(X = counts, fit = tm3_r1_init)

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
tm4_r1_init <- fit_poisson_nmf(
  X = counts, 
  fit0 = tm4_fit0, 
  numiter = 1000,
  control = list(nc = 7)
  )

ll_tm4 <- sum(loglik_poisson_nmf(X = counts, fit = tm4_r1_init))
ll_per_pt_tm4 <- loglik_poisson_nmf(X = counts, fit = tm4_r1_init)

############# Rank 3 log1p model with Rank 1 Initialization ###########
library(log1pNMF)

set.seed(10)
log1p_fit3 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 3,
  loglik = "exact",
  control = list(maxiter = 1000)
)

ll_log1p_fit3 <- logLik(log1p_fit3, Y = counts)

ll_per_pt_log1p_fit3 <- rowSums(dpois(
  x = counts,
  lambda = fitted(log1p_fit3),
  log = TRUE
))

set.seed(10)
log1p_fit4 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 4,
  loglik = "exact",
  control = list(maxiter = 1000)
)

ll_log1p_fit4 <- logLik(log1p_fit4, Y = counts)

ll_per_pt_log1p_fit4 <- rowSums(dpois(
  x = counts,
  lambda = fitted(log1p_fit4),
  log = TRUE
))

ll_df <- data.frame(
  group = rep(samples$label, 4),
  loglik = c(
    ll_per_pt_log1p_fit3,
    ll_per_pt_log1p_fit4,
    ll_per_pt_tm3,
    ll_per_pt_tm4
  ),
  rank = c(
    rep(3, nrow(samples)),
    rep(4, nrow(samples)),
    rep(3, nrow(samples)),
    rep(4, nrow(samples))
  ),
  model = c(
    rep("log1p", 2 * nrow(samples)),
    rep("topic model", 2 * nrow(samples))
  )
)

ggplot(ll_df, aes(x = group, y = loglik)) +
  geom_boxplot() +
  facet_grid(rows  = vars(rank),
             cols  = vars(model)) +
  labs(
    x     = "Group",
    y     = "Log-likelihood"
  ) +
  theme_cowplot()


# want to look at structure plots to understand what's really going on here
topic_colors <- c("olivedrab", "tomato","darkblue","dodgerblue")
n <- nrow(samples)
L0 <- poisson2multinom(tm4_r1_init)$L
L <- L0
L0[,2] <- L[,4]
L0[,4] <- L[,2]
g2 <- structure_plot(L0,grouping = samples$label,topics = c(2,3,4,1),
                     loadings_order = 1:n,colors = topic_colors)$plot +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5)) + 
  ggtitle("Topic Model") 


