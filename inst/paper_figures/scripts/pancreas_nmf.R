library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)
library(log1pNMF)

set.seed(1)

load("~/Downloads/pancreas.RData")

i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]
sample_info <- transform(sample_info,celltype = factor(celltype))

genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]

n <- nrow(counts)
p <- ncol(counts)

K <- 9

library(fastTopics)
nmf_k1 <- fastTopics:::fit_pnmf_rank1(X = counts)
L0 <- nmf_k1$L %>%
  cbind(
    matrix(
      data = 1e-10,
      nrow = n,
      ncol = K - 1
    )
  )
F0 <- nmf_k1$F %>%
  cbind(
    matrix(
      data = 1e-10,
      nrow = p,
      ncol = K - 1
    )
  )

rownames(L0) <- rownames(counts)
rownames(F0) <- colnames(counts)

nmf_fit0 <- fastTopics::init_poisson_nmf(
  X = counts,
  F = F0,
  L = L0
)
nmf <- fit_poisson_nmf(X = counts, fit0 = nmf_fit0, control = list(nc = 7), numiter = 250)

readr::write_rds(nmf, "~/Documents/data/log1pNMF/pancreas/pancreas_pois_nmf_k9_exact_250_iter.rds")

celltype <- sample_info$celltype
celltype <-
  factor(celltype,
         c("acinar","ductal","activated_stellate","quiescent_stellate",
           "endothelial","macrophage","mast","schwann","alpha","beta",
           "delta","gamma","epsilon"))

levels(celltype)[levels(celltype) == "activated_stellate"] <- "activated PSC"
levels(celltype)[levels(celltype) == "quiescent_stellate"] <- "quiescent PSC"

structure_plot(nmf, grouping = celltype,gap = 30,perplexity = 70,n = Inf, font.size = 12)

