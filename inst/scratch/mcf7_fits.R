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

set.seed(1)
fgpca_fit <- fit_glmpca_pois(
  Y = t(counts),
  K = 2,
  control = list(maxiter = 1000)
)

fit_list[["glmpca"]] <- fgpca_fit
readr::write_rds(fit_list, "~/Documents/data/mcf7_fit_list.rds")

