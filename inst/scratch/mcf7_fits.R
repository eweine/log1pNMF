library(rsvd)
library(fastglmpca)
library(fastTopics)
library(flashier)
library(data.table)
library(GEOquery)
library(ggplot2)
library(cowplot)
library(log1pNMF)

genes <- fread("../data/Human.GRCh38.p13.annot.tsv.gz",
               sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(genes) <- "data.frame"
genes <- genes[1:10]
genes <- transform(genes,
                   GeneType = factor(GeneType),
                   Status   = factor(Status))

counts <- fread("../data/GSE152749_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(counts) <- "data.frame"
rownames(counts) <- counts$GeneID
counts <- counts[,-1]
counts <- as.matrix(counts)
storage.mode(counts) <- "double"
counts <- t(counts)
ids <- rownames(counts)

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

K <- 3
cc_vec <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000)

fit_list <- list()

for (cc in cc_vec) {

  set.seed(10)
  fit <- fit_poisson_log1p_nmf(
    Y = counts,
    K = K,
    cc = cc,
    loglik = "exact",
    control = list(maxiter = 1000)
  )

  fit_list[[as.character(cc)]] <- fit

}

# finally, fit a topic model with a rank 1 initialization here
r1_fit <- fastTopics:::fit_pnmf_rank1(counts)

init_LL <- cbind(
  r1_fit$L,
  matrix(
    data = 1e-8,
    nrow = nrow(counts),
    ncol = K - 1
  )
)
rownames(init_LL) <- rownames(counts)

init_FF <- cbind(
  r1_fit$F,
  matrix(
    data = 1e-8,
    nrow = ncol(counts),
    ncol = K - 1
  )
)
rownames(init_FF) <- colnames(counts)

fit0 <- init_poisson_nmf(X = counts, L = init_LL, F = init_FF)

nmf_fit <- fit_poisson_nmf(
  X = counts,
  fit0 = fit0,
  numiter = 1000,
  control = list(nc = 7)
)

fit_list[[as.character(Inf)]] <- nmf_fit

readr::write_rds(fit_list, "~/Documents/data/mcf7_fit_list.rds")
fit_list <- readr::read_rds("~/Documents/data/mcf7_fit_list.rds")

hoyer_sparsity <- function(x) {

  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))

}

for (cc in cc_vec) {

  fit <- fit_list[[as.character(cc)]]
  fit_list[[as.character(cc)]]$l_sparsity <- median(apply(
    fit$LL, 2, hoyer_sparsity
  ))
  fit_list[[as.character(cc)]]$f_sparsity <- median(apply(
    fit$FF, 2, hoyer_sparsity
  ))

  fit_list[[as.character(cc)]]$f_cor <- median(
    abs(cor(fit$FF, method = "spearman"))[lower.tri(diag(K))]
  )
  
  fit_list[[as.character(cc)]]$ll <- logLik(
    fit, Y = counts
  )

}

fit <- fit_list[[as.character(Inf)]]

fit_list[[as.character(Inf)]]$ll <- sum(loglik_poisson_nmf(X = counts, fit))

fit_list[[as.character(Inf)]]$l_sparsity <- median(apply(
  fit$L, 2, hoyer_sparsity
))
fit_list[[as.character(Inf)]]$f_sparsity <- median(apply(
  fit$F, 2, hoyer_sparsity
))


fit_list[[as.character(Inf)]]$f_cor <- median(
  abs(cor(fit$F, method = "spearman"))[lower.tri(diag(K))]
)

l_sparsity_vec <- unlist(lapply(fit_list, function(x) {x$l_sparsity}))
f_sparsity_vec <- unlist(lapply(fit_list, function(x) {x$f_sparsity}))
f_cor_vec <- unlist(lapply(fit_list, function(x) {x$f_cor}))
ll_vec <- unlist(lapply(fit_list, function(x) {x$ll}))

library(dplyr)
df_sparsity_l <- data.frame(
  cc = as.numeric(names(l_sparsity_vec)),
  sparsity = l_sparsity_vec
) %>% filter(is.finite(cc))


df_ll <- data.frame(
  cc = as.numeric(names(ll_vec)),
  ll = ll_vec
) %>% filter(is.finite(cc))

df_sparsity_f <- data.frame(
  cc = as.numeric(names(f_sparsity_vec)),
  sparsity = f_sparsity_vec
) %>% filter(is.finite(cc))

ggplot(data = df_ll, aes(x = cc, y = ll)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  xlab("c (log10 scale)") +
  ylab("loglik") +
  geom_hline(yintercept = ll_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.004, y=ll_vec["Inf"] + 10000, label="ID Link", color="red"
  )


ggplot(data = df_sparsity_l, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  xlab("c (log10 scale)") +
  ylab("Median Loading Sparsity") +
  geom_hline(yintercept = l_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.004, y=l_sparsity_vec["Inf"] + 0.05, label="ID Link", color="red"
  )

ggplot(data = df_sparsity_f, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  xlab("c (log10 scale)") +
  ylab("Median Factor Sparsity") +
  geom_hline(yintercept = f_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.004, y=f_sparsity_vec["Inf"] + 0.05, label="ID Link", color="red"
  )

df_cor <- data.frame(
  cc = as.numeric(names(f_cor_vec)),
  correlation = f_cor_vec
) %>% filter(is.finite(cc))

ggplot(data = df_cor, aes(x = cc, y = correlation)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  xlab("c (log10 scale)") +
  ylab("Median Abs. Factor Correlation") +
  geom_hline(yintercept = f_cor_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.004, y=f_cor_vec["Inf"] + 0.05, label="ID Link", color="red"
  )

library(log1pNMF)
topic_colors <- c("olivedrab","dodgerblue","darkblue","tomato")
normalized_structure_plot(fit_list$`1e-04`,grouping = samples$label,topics = 4:1,
               colors = topic_colors)$plot +
  labs(y = "membership")

structure_plot(fit, grouping = samples$label,topics = 4:1,
               colors = topic_colors)$plot + labs(y = "membership")

