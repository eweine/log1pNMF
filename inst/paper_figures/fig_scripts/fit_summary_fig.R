# here, want to make a figure encompassing all datasets comparing sparsity
# and orthogonality

load("../data/experiment_results.Rdata")
dat <- readr::read_csv("../data/raw_data/bbc_news_text_complexity_summarization.csv")

set.seed(1)
library(fastTopics)
library(tm)
library(SnowballC)

my_corpus <- VCorpus(VectorSource(dat$text))

addspace <- content_transformer(function(x, pattern) {
  return(gsub(pattern, " ", x))
})

my_corpus <- tm_map(my_corpus, addspace, "-")

my_corpus <- tm_map(my_corpus, removeNumbers)

# Transform to lower case (need to wrap in content_transformer)
my_corpus <- tm_map(my_corpus,content_transformer(tolower))
my_corpus <- tm_map(my_corpus, removeWords, stopwords("SMART"))
my_corpus <- tm_map(my_corpus, removePunctuation)
my_corpus <- tm_map(my_corpus, stripWhitespace)

my_corpus <- tm_map(my_corpus, stemDocument)

dtm <- DocumentTermMatrix(my_corpus)
dtm2 <- Matrix::sparseMatrix(
  i = dtm$i,
  j = dtm$j,
  x = dtm$v
)


colnames(dtm2) <- dtm$dimnames$Terms
words_to_use <- which(Matrix::colSums(dtm2>0)>4)

dtm2 <- dtm2[,words_to_use]

s <- Matrix::rowSums(dtm2)
s <- s / mean(s)

library(log1pNMF)
library(Matrix)
library(dplyr)
#cc_vec <- c(1e3)

n <- nrow(dtm2)
p <- ncol(dtm2)

K <- 10

hoyer_sparsity <- function(x) {
  
  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))
  
}

fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

for (cc in cc_vec) {
  
  fit_list[[as.character(cc)]] <- res_list$bbc[[as.character(cc)]]
  
  colnames(fit_list[[as.character(cc)]]$LL) <- paste0("k", 1:K)
  
  fit_list[[as.character(cc)]]$l_sparsity <- apply(
    fit_list[[as.character(cc)]]$LL, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$f_sparsity <- apply(
    fit_list[[as.character(cc)]]$FF, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$final_loglik <- logLik(
    object = fit_list[[as.character(cc)]], Y = dtm2
  )
  
  fit_list[[as.character(cc)]]$cor_mat <- cor(fit_list[[as.character(cc)]]$FF, method = "spearman")
  
}

fit_list[["Inf"]] <- res_list$bbc$`Inf`
fit_list[["Inf"]]$final_loglik <- NA

fit_list[["Inf"]]$Ls <- Matrix::Diagonal(x = 1/s) %*% fit_list[["Inf"]]$L

fit_list[["Inf"]]$l_sparsity <- apply(
  fit_list[["Inf"]]$Ls, 2, hoyer_sparsity
)

fit_list[["Inf"]]$f_sparsity <- apply(
  fit_list[["Inf"]]$F, 2, hoyer_sparsity
)

fit_list[["Inf"]]$cor_mat <- cor(fit_list[["Inf"]]$F, method = "spearman")

l_sparsity_vec <- unlist(lapply(fit_list, function(x) {mean(x$l_sparsity)}))
f_sparsity_vec <- unlist(lapply(fit_list, function(x) {mean(x$f_sparsity)}))
cor_vec <- unlist(
  lapply(fit_list, function(x) {mean(abs(x$cor_mat[lower.tri(x$cor_mat)]))})
)

loglik_vec <- unlist(lapply(fit_list, function(x) {x$final_loglik}))

df_sparsity_l <- data.frame(
  cc = as.numeric(names(l_sparsity_vec)),
  sparsity = l_sparsity_vec
) %>% filter(is.finite(cc))

df_sparsity_f <- data.frame(
  cc = as.numeric(names(f_sparsity_vec)),
  sparsity = f_sparsity_vec
) %>% filter(is.finite(cc))

df_cor <- data.frame(
  cc = as.numeric(names(cor_vec)),
  correlation = cor_vec
) %>% filter(is.finite(cc))

library(ggpubr)
library(ggplot2)

g1 <- ggplot(data = df_cor, aes(x = cc, y = correlation)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c") +
  ylab("Mean Abs. Factor Correlation") +
  geom_hline(yintercept = cor_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.008, y=cor_vec["Inf"] + 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

g2 <- ggplot(data = df_sparsity_l, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c") +
  ylab("Mean Loading Sparsity") +
  geom_hline(yintercept = l_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.008, y=l_sparsity_vec["Inf"] + 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

g3 <- ggplot(data = df_sparsity_f, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c") +
  ylab("Mean Factor Sparsity") +
  geom_hline(yintercept = f_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.008, y=f_sparsity_vec["Inf"] + 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

rm(list = setdiff(ls(), c("g1", "g2", "g3", "hoyer_sparsity", "res_list")))

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

K <- 3
x <- colSums(counts > 0)
i <- which(x > 3 &
             genes$GeneType == "protein-coding" &
             genes$Status == "active")
genes  <- genes[i,]
counts <- counts[,i]

s <- Matrix::rowSums(counts)
s <- s / mean(s)

fit_list <- list()
cc_vec <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4)

for (cc in cc_vec) {
  
  fit_list[[as.character(cc)]] <- res_list$mcf7[[as.character(cc)]]
  
  colnames(fit_list[[as.character(cc)]]$LL) <- paste0("k", 1:K)
  
  fit_list[[as.character(cc)]]$l_sparsity <- apply(
    fit_list[[as.character(cc)]]$LL, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$f_sparsity <- apply(
    fit_list[[as.character(cc)]]$FF, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$final_loglik <- logLik(
    object = fit_list[[as.character(cc)]], Y = counts
  )
  
  fit_list[[as.character(cc)]]$cor_mat <- cor(fit_list[[as.character(cc)]]$FF, method = "spearman")
  
}

fit_list[["Inf"]] <- res_list$mcf7$`Inf`
fit_list[["Inf"]]$final_loglik <- NA

fit_list[["Inf"]]$Ls <- Matrix::Diagonal(x = 1/s) %*% fit_list[["Inf"]]$L

fit_list[["Inf"]]$l_sparsity <- apply(
  fit_list[["Inf"]]$Ls, 2, hoyer_sparsity
)

fit_list[["Inf"]]$f_sparsity <- apply(
  fit_list[["Inf"]]$F, 2, hoyer_sparsity
)

fit_list[["Inf"]]$cor_mat <- cor(fit_list[["Inf"]]$F, method = "spearman")

l_sparsity_vec <- unlist(lapply(fit_list, function(x) {mean(x$l_sparsity)}))
f_sparsity_vec <- unlist(lapply(fit_list, function(x) {mean(x$f_sparsity)}))
cor_vec <- unlist(
  lapply(fit_list, function(x) {mean(abs(x$cor_mat[lower.tri(x$cor_mat)]))})
)

loglik_vec <- unlist(lapply(fit_list, function(x) {x$final_loglik}))

df_sparsity_l <- data.frame(
  cc = as.numeric(names(l_sparsity_vec)),
  sparsity = l_sparsity_vec
) %>% filter(is.finite(cc))

df_sparsity_f <- data.frame(
  cc = as.numeric(names(f_sparsity_vec)),
  sparsity = f_sparsity_vec
) %>% filter(is.finite(cc))

df_cor <- data.frame(
  cc = as.numeric(names(cor_vec)),
  correlation = cor_vec
) %>% filter(is.finite(cc))

library(ggpubr)
library(ggplot2)

g4 <- ggplot(data = df_cor, aes(x = cc, y = correlation)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  xlab("c") +
  ylab("Mean Abs. Factor Correlation") +
  geom_hline(yintercept = cor_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.003, y=cor_vec["Inf"] + 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

g5 <- ggplot(data = df_sparsity_l, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  xlab("c") +
  ylab("Mean Loading Sparsity") +
  geom_hline(yintercept = l_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.003, y=l_sparsity_vec["Inf"] + 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

g6 <- ggplot(data = df_sparsity_f, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  xlab("c") +
  ylab("Mean Factor Sparsity") +
  geom_hline(yintercept = f_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.003, y=f_sparsity_vec["Inf"] + 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

rm(list = setdiff(ls(), c(
  "g1", "g2", "g3", "g4", "g5", "g6", "hoyer_sparsity", "res_list"
  )))

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")
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

s <- Matrix::rowSums(counts)
s <- s / mean(s)

K <- 13

fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

for (cc in cc_vec) {
  
  fit_list[[as.character(cc)]] <- res_list$pancreas[[as.character(cc)]]
  
  colnames(fit_list[[as.character(cc)]]$LL) <- paste0("k", 1:K)
  
  fit_list[[as.character(cc)]]$l_sparsity <- apply(
    fit_list[[as.character(cc)]]$LL, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$f_sparsity <- apply(
    fit_list[[as.character(cc)]]$FF, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$final_loglik <- logLik(
    object = fit_list[[as.character(cc)]], Y = counts
  )
  
  fit_list[[as.character(cc)]]$cor_mat <- cor(fit_list[[as.character(cc)]]$FF, method = "spearman")
  
}

fit_list[["Inf"]] <- res_list$pancreas$`Inf`
fit_list[["Inf"]]$final_loglik <- NA

fit_list[["Inf"]]$Ls <- Matrix::Diagonal(x = 1/s) %*% fit_list[["Inf"]]$L

fit_list[["Inf"]]$l_sparsity <- apply(
  fit_list[["Inf"]]$Ls, 2, hoyer_sparsity
)

fit_list[["Inf"]]$f_sparsity <- apply(
  fit_list[["Inf"]]$F, 2, hoyer_sparsity
)

fit_list[["Inf"]]$cor_mat <- cor(fit_list[["Inf"]]$F, method = "spearman")

l_sparsity_vec <- unlist(lapply(fit_list, function(x) {mean(x$l_sparsity)}))
f_sparsity_vec <- unlist(lapply(fit_list, function(x) {mean(x$f_sparsity)}))
cor_vec <- unlist(
  lapply(fit_list, function(x) {mean(abs(x$cor_mat[lower.tri(x$cor_mat)]))})
)

loglik_vec <- unlist(lapply(fit_list, function(x) {x$final_loglik}))

df_sparsity_l <- data.frame(
  cc = as.numeric(names(l_sparsity_vec)),
  sparsity = l_sparsity_vec
) %>% filter(is.finite(cc))

df_sparsity_f <- data.frame(
  cc = as.numeric(names(f_sparsity_vec)),
  sparsity = f_sparsity_vec
) %>% filter(is.finite(cc))

df_cor <- data.frame(
  cc = as.numeric(names(cor_vec)),
  correlation = cor_vec
) %>% filter(is.finite(cc))

g7 <- ggplot(data = df_cor, aes(x = cc, y = correlation)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c") +
  ylab("Mean Abs. Factor Correlation") +
  geom_hline(yintercept = cor_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.008, y=cor_vec["Inf"] + 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

g8 <- ggplot(data = df_sparsity_l, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c") +
  ylab("Mean Loading Sparsity") +
  geom_hline(yintercept = l_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.008, y=l_sparsity_vec["Inf"] + 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

g9 <- ggplot(data = df_sparsity_f, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c") +
  ylab("Mean Factor Sparsity") +
  geom_hline(yintercept = f_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.008, y=f_sparsity_vec["Inf"] - 0.05, label="c = \u221E", color="red",
    size = 5
  ) +
  theme(
    axis.title.x = element_text(size = 14.5),
    axis.title.y = element_text(size = 14.5)
  )

mcf7_sum_fig <- ggarrange(g4,g5,g6, nrow = 1, labels = c("A", "B", "C"))
mcf7_sum_fig <- annotate_figure(mcf7_sum_fig,
                     top = text_grob("MCF-7 K = 3", size = 25, face = "bold"))

bbc_sum_fig <- ggarrange(g1,g2,g3, nrow = 1, labels = c("G", "H", "I"))
bbc_sum_fig <- annotate_figure(bbc_sum_fig,
                                top = text_grob("BBC K = 10", size = 25, face = "bold"))


lsa_sum_fig <- ggarrange(g7,g8,g9, nrow = 1, labels = c("D", "E", "F"))
lsa_sum_fig <- annotate_figure(lsa_sum_fig,
                               top = text_grob("Pancreas K = 13", size = 25, face = "bold"))

all_dataset_sum_fig <- ggarrange(mcf7_sum_fig, lsa_sum_fig, bbc_sum_fig, nrow = 3, ncol = 1)

ggsave(
  "../images/fit_summary.png",
  all_dataset_sum_fig,
  device = "png",
  width = 12.66,
  height = 12.66
)
