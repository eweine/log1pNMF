library(log1pNMF)
library(Matrix)
library(dplyr)
set.seed(1)

counts <- readr::read_rds("heart.rds")
reticulate::use_condaenv("nmf")
ad <- reticulate::import("anndata")
dat <- ad$read_h5ad("/Users/eweine/Downloads/healthy_human_4chamber_map_unnormalized_V3.h5ad")
obs <- dat$obs
obs <- obs %>%
  dplyr::filter(rownames(obs) %in% rownames(counts))

obs <- obs %>%
  dplyr::filter(!(Cluster %in% c("08. Macrophage", "16. Neuronal", "17. Lymphocyte")))

counts <- counts[rownames(obs), ]

tot <- obs %>% dplyr::group_by(Cluster) %>%
  dplyr::summarise(tot = dplyr::n())

s <- Matrix::rowSums(counts)
s <- s / mean(s)

genes_to_use <- which(Matrix::colSums(counts>0)>9)
counts <- counts[,genes_to_use]

counts <- counts[, !grepl("^MT-|^RPS|^RPL|^MALAT1$", colnames(counts))]

K <- 12
cc_vec <- c(1e-3, 1)
n_iter <- 100

for (cc in cc_vec) {
  
  print(cc)
  
  set.seed(1)
  fit <- fit_poisson_log1p_nmf(
    Y = counts,
    K = K,
    s = s,
    cc = cc,
    loglik = "exact",
    init_method = "rank1",
    control = list(
      maxiter = n_iter,
      threads = 35
    )
  )
  
  readr::write_rds(
    fit, glue::glue("heart_log1p_c{cc}_k{K}_exact_{n_iter}_iter.rds")
  )
  
}

n <- nrow(counts)
p <- ncol(counts)

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
nmf <- fit_poisson_nmf(X = counts, fit0 = nmf_fit0, control = list(nc = 7), numiter = n_iter)

reticulate::use_condaenv("nmf")

ad <- reticulate::import("anndata")
structure_plot(nmf, grouping = obs$Cluster, gap = 20)


