high_expressed_genes <- 500
low_expressed_genes <- 500

p <- high_expressed_genes + low_expressed_genes
n <- 1000

LL <- matrix(
  data = c(
    rep(1, 333), rep(0, 1000 - 333),
    rep(0, 333), rep(1, 333), rep(0, 334),
    rep(0, 1000 - 334), rep(1, 334)
  ),
  nrow = n,
  ncol = 3
)


FF <- matrix(
  data = c(
    rep(3, round(low_expressed_genes / 2)), 
    rep(0, round(low_expressed_genes / 2)),
    rep(1100, round(high_expressed_genes / 2)),
    rep(1000, round(high_expressed_genes / 2)),
    rep(0, round(low_expressed_genes / 2)), 
    rep(3, round(low_expressed_genes / 2)),
    rep(1000, round(high_expressed_genes / 2)),
    rep(1100, round(high_expressed_genes / 2)),
    rep(0, round(low_expressed_genes / 2)), 
    rep(0, round(low_expressed_genes / 2)),
    rep(1051.5, round(high_expressed_genes / 2)),
    rep(1051.5, round(high_expressed_genes / 2))
  ),
  nrow = p,
  ncol = 3
)

Lambda <- tcrossprod(LL, FF)
set.seed(1)
Y <- matrix(
  data = rpois(n = n * p, lambda = as.vector(Lambda)),
  nrow = n,
  ncol = p
)

library(fastTopics)
rownames(Y) <- paste0("cell", 1:n)
colnames(Y) <- paste0("gene", 1:p)

ft_r1 <- fastTopics:::fit_pnmf_rank1(Y)
init_LL <- cbind(
  ft_r1$L,
  matrix(
    data = 1e-3,
    nrow = n,
    ncol = 2
  )
)
rownames(init_LL) <- rownames(Y)

init_FF <- cbind(
  ft_r1$F,
  matrix(
    data = 1e-3,
    nrow = p,
    ncol = 2
  )
)
rownames(init_FF) <- colnames(Y)

ft_init <- init_poisson_nmf(X = Y, F = init_FF, L = init_LL)

ft_mod <- fit_poisson_nmf(X = Y, fit0 = ft_init)

group <- c(rep("A", 333), rep("B", 333), rep("C", 334))


structure_plot(log1pNMF:::normalize_bars(ft_mod$L), grouping = group, gap = 10, loadings_order = 1:n)

library(log1pNMF)
set.seed(1)
log1p_mod <- fit_poisson_log1p_nmf(
  Y = Y, K = 3, loglik = "exact",
  control = list(maxiter = 250)
)

normalized_structure_plot(
  log1p_mod, 
  grouping = group, 
  gap = 10,
  loadings_order = 1:n
  )

library(fastglmpca)
set.seed(1)
gpca_fit <- fit_glmpca_pois(
  Y = t(Y),
  K = 2,
  control = list(maxiter = 1000)
)

gpca_df <- as.data.frame(gpca_fit$V)
gpca_df$group <- group

ggplot(data = gpca_df, aes(x = k_1, y = k_2)) +
  geom_point(aes(color = group))

