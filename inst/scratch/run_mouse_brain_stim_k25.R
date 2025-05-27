load("~/Downloads/mouse_brain_stim.Rdata")

library(dplyr)

set.seed(1)
cells <- cells %>% dplyr::filter(!is.na(maintype))
counts <- counts[rownames(counts) %in% cells$...1, ]
counts <- counts[Matrix::rowSums(counts) > 0, ]

genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[, genes_to_use]
counts <- as(counts, "CsparseMatrix")

n <- nrow(counts)
p <- ncol(counts)
K <- 15

library(fastTopics)
pois_nmf_k1 <- fastTopics:::fit_pnmf_rank1(counts)

set.seed(1)
init_LL <- pois_nmf_k1$L %>%
  cbind(
    matrix(
      data = 1e-10,
      nrow = n,
      ncol = K - 1
    )
  )

init_FF <- pois_nmf_k1$F %>%
  cbind(
    matrix(
      data = 1e-10,
      nrow = p,
      ncol = K - 1
    )
  )

rownames(init_LL) <- rownames(counts)
rownames(init_FF) <- colnames(counts)

pois_nmf_k15 <- init_poisson_nmf(
  X = counts,
  L = init_LL,
  F = init_FF
)

pois_nmf_k15 <- fit_poisson_nmf(
  X = counts,
  fit0 = pois_nmf_k15,
  control = list(nc = 7),
  numiter = 250
)

readr::write_rds(pois_nmf_k15, "~/Documents/data/log1pNMF/experiment_results/mouse_brain_k12_nmf_pois.rds")

# the last thing to do here is to run frobenius NMF on the log1p
# transformed data. I think it is best to do this on Midway
# once I get home

