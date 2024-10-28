load("~/Downloads/mouse_brain_stim.Rdata")

library(dplyr)

set.seed(1)
cells <- cells %>% dplyr::filter(!is.na(maintype))
counts <- counts[rownames(counts) %in% cells$...1, ]
counts <- counts[Matrix::rowSums(counts) > 0, ]

counts <- counts[, Matrix::colSums(counts) > 0]
counts <- as(counts, "CsparseMatrix")

n <- nrow(counts)
p <- ncol(counts)
K <- 25

rs <- Matrix::rowSums(counts)
s <- rs / mean(rs)

set.seed(1)
log1p_k1 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 1,
  maxiter = 10,
  approx_range = c(0, 1.25),
  s = s
)

init_LL <- log1p_k1$U %>%
  cbind(
    matrix(
      data = rexp(
        n = n * (K - 1), rate = 15
      ),
      nrow = n,
      ncol = K - 1
    )
  )

init_FF <- log1p_k1$V %>%
  cbind(
    matrix(
      data = rexp(
        n = p * (K - 1), rate = 15
      ),
      nrow = p,
      ncol = K - 1
    )
  )

set.seed(1)
log1p_k25 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = K,
  maxiter = 100,
  approx_range = c(0, 1.25),
  s = s,
  init_U = init_LL,
  init_V = init_FF
)

readr::write_rds(log1p_k25, "~/Documents/data/passPCA/experiment_results/mouse_brain_k25_log1p_pois.rds")

library(fastTopics)
pois_nmf_k1 <- fastTopics:::fit_pnmf_rank1(counts)

set.seed(1)
init_LL <- pois_nmf_k1$L %>%
  cbind(
    matrix(
      data = rexp(
        n = n * (K - 1), rate = 15
      ),
      nrow = n,
      ncol = K - 1
    )
  )

init_FF <- pois_nmf_k1$F %>%
  cbind(
    matrix(
      data = rexp(
        n = p * (K - 1), rate = 15
      ),
      nrow = p,
      ncol = K - 1
    )
  )

rownames(init_LL) <- rownames(counts)
rownames(init_FF) <- colnames(counts)

pois_nmf_k25 <- init_poisson_nmf(
  X = counts,
  L = init_LL,
  F = init_FF
)

pois_nmf_k25 <- fit_poisson_nmf(
  X = counts,
  fit0 = pois_nmf_k25,
  control = list(nc = 7)
)

readr::write_rds(pois_nmf_k25, "~/Documents/data/passPCA/experiment_results/mouse_brain_k25_nmf_pois.rds")

# the last thing to do here is to run frobenius NMF on the log1p
# transformed data. I think it is best to do this on Midway
# once I get home

