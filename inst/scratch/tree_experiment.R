# here, I'm interested in generating some sort of hierarchical data
# I think the simplest form of hierarchy is a multi-group model, where
# there are n groups. All groups are loaded on a "baseline" factor
# and then each group is loaded on a sparse change factor.

# number of genes
p <- 500
# number of cells
n <- 1000
K <- 7
LL <- matrix(
  data = 0,
  nrow = n,
  ncol = K
)

FF <- matrix(
  data = 0,
  nrow = p,
  ncol = K
)

set.seed(1)
# everyone loaded on the first factor
LL[,1] <- 1
LL[1:500, 2] <- 1
LL[501:1000, 3] <- 1
LL[1:250, 4] <- 1
LL[251:500, 5] <- 1
LL[501:750, 6] <- 1
LL[751:1000, 7] <- 1

# base rate
FF[,1] <- 0.5
all_idx <- sample(1:500)  # Shuffle all indices once
split_idx <- split(all_idx, rep(1:K, each = 25))  # Divide into K groups

for (k in 2:K) {
  FF[split_idx[[k]], k] <- FF[split_idx[[k]], k] + 3
}

Lambda <- tcrossprod(LL, FF)

set.seed(1)
Y <- matrix(
  data = rpois(
    n = n * p,
    lambda = as.vector(Lambda)
  ),
  nrow = n,
  ncol = p
)
Y <- as(Y, "CsparseMatrix")

set.seed(1)
log1p_ft7 <- fit_poisson_log1p_nmf(
  Y = Y, K = K, s = FALSE, cc = 1e-3, loglik = "exact"
)

normalized_structure_plot(log1p_ft7, loadings_order = 1:n)
library(fastTopics)
structure_plot(LL, loadings_order = 1:n)
normalized_structure_plot(log1p_ft7, loadings_order = 1:n)


library(fastTopics)
set.seed(1)
ft4 <- fit_poisson_nmf(X = Y, k = 4, control = list(nc = 7), numiter = 250)

set.seed(1)
ft7 <- fit_poisson_nmf(X = Y, k = 7, control = list(nc = 7), numiter = 250)

library(log1pNMF)

structure_plot(ft4, loadings_order = 1:n)
structure_plot(ft7, loadings_order = 1:n)

nmf_k1 <- fastTopics:::fit_pnmf_rank1(Y)
library(dplyr)
nmf_LL <- nmf_k1$L %>%
  cbind(
    matrix(
      data = rexp(
        n = n * (K - 1), rate = 15
      ),
      nrow = n,
      ncol = K - 1
    )
  )
rownames(nmf_LL) <- rownames(Y)

nmf_FF <- nmf_k1$F %>%
  cbind(
    matrix(
      data = rexp(
        n = p * (K - 1), rate = 15
      ),
      nrow = p,
      ncol = K - 1
    )
  )

rownames(nmf_FF) <- colnames(Y)

set.seed(1)
nmf_fit0 <- fastTopics::init_poisson_nmf(
  X = Y,
  L = nmf_LL,
  F = nmf_FF
)

nmf_fit <- fastTopics::fit_poisson_nmf(
  X = Y,
  fit0 = nmf_fit0,
  control = list(nc = 7), numiter = 250
)

structure_plot(nmf_fit, loadings_order = 1:n)
