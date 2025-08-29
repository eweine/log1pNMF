n_cells <- 1000
n_genes <- 1000

library(fastTopics)
library(dplyr)

set.seed(1)
LL <- matrix(
  data = c(
    rep(1, 200), rep(0, 800),
    rep(0, 200), rep(1, 200), rep(0, 600),
    rep(0, 400), rep(1, 200), rep(0, 400),
    rep(0, 600), rep(1, 200), rep(0, 200),
    rep(0, 800), rep(1, 200),
    rbinom(n = 1000, size = 1, prob = 0.5)
  ),
  nrow = n_cells,
  ncol = 6
)

LL_df <- data.frame(
  idx = 1:1000,
  k1 = LL[,1],
  k2 = LL[,2],
  k3 = LL[,3],
  k4 = LL[,4],
  k5 = LL[,5],
  k6 = LL[,6]
) %>% 
  dplyr::arrange(
    desc(k1), desc(k2), desc(k3), desc(k4), desc(k5), desc(k6)
  )

LL <- LL[LL_df$idx, ]

structure_plot(LL, loadings_order = 1:n_cells, topics = rev(1:6))

mult_vec <- rep(2.5, 1000)
mult_vec[sample(1:1000, size = 800)] <- 0

FF <- matrix(
  data = c(
    rexp(200), rep(0, 800),
    rep(0, 200), rexp(200), rep(0, 600),
    rep(0, 400), rexp(200), rep(0, 400),
    rep(0, 600), rexp(200), rep(0, 200),
    rep(0, 800), rexp(200),
    mult_vec
  ),
  nrow = n_genes,
  ncol = 6
)

Lambda <- expm1(tcrossprod(LL, FF))

Y <- matrix(
  data = rpois(n = n_cells * n_genes, lambda = as.vector(Lambda)),
  nrow = n_cells,
  ncol = n_genes
)
Y <- as(Y, "CsparseMatrix")
Y <- Y[,Matrix::colSums(Y) > 0]

library(fastTopics)

set.seed(1)
ft_fit <- fit_poisson_nmf(X = Y, k = 6)
structure_plot(ft_fit, loadings_order = 1:n_cells)

library(log1pNMF)

set.seed(1)
log1p_fit <- fit_poisson_log1p_nmf(
  Y = Y, K = 6, loglik = "exact", init_method = "random",
  control = list(maxiter = 250)
)
normalized_structure_plot(log1p_fit, loadings_order = 1:n_cells)




