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

K <- 4

set.seed(1)
fit0 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e4,
                              loglik = "exact",init_method = "rank1",
                              control = list(maxiter = 20,verbose = TRUE))

set.seed(1)
fit1 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e4,
                              loglik = "exact",init_method = "rank1",
                              control = list(maxiter = 10,verbose = TRUE))
print(any(is.na(fit1$objective_trace)))

fit2 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e4,loglik = "exact",
                              init_LL = fit1$LL,init_FF = fit1$FF,
                              control = list(maxiter = 10,verbose = TRUE))
print(any(is.na(fit2$objective_trace)))
