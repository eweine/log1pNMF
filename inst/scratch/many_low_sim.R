n_genes <- 1000
n_groups <- 4
n_cells_per_group <- 250

high_expr <- 50
low_expr <- 1.5

lambda_list <- list()
Y_list <- list()

lambda_list[[1]] <- c(
  rep(low_expr, 960),
  rep(high_expr, 10),
  rep(0, 30)
)

lambda_list[[2]] <- c(
  rep(low_expr, 960),
  rep(0, 10),
  rep(high_expr, 10),
  rep(0, 20)
)

lambda_list[[3]] <- c(
  rep(low_expr, 960),
  rep(0, 20),
  rep(high_expr, 10),
  rep(0, 10)
)

lambda_list[[4]] <- c(
  rep(low_expr, 960),
  rep(0, 30),
  rep(high_expr, 10)
)

set.seed(1)

for (group in 1:n_groups) {
  
  Lambda <- matrix(
    data = rep(lambda_list[[group]], n_cells_per_group),
    nrow = n_cells_per_group,
    ncol = n_genes,
    byrow = TRUE
  )
  
  Y <- matrix(
    data = rpois(n = n_cells_per_group * n_genes, lambda = as.vector(Lambda)),
    nrow = n_cells_per_group,
    ncol = n_genes
  )
  
  Y_list[[group]] <- Y
  
}

Y <- do.call(rbind, Y_list)

library(log1pNMF)
set.seed(1)
log1p_mod <- fit_poisson_log1p_nmf(
  Y = Y,
  K = 4,
  loglik = "exact",
  init_method = "random",
  control = list(maxiter = 1000),
  cc = 1e-3
)

normalized_structure_plot(log1p_mod, loadings_order = 1:1000)

set.seed(1)
log1p_mod_r1_init <- fit_poisson_log1p_nmf(
  Y = Y,
  K = 4,
  loglik = "exact",
  init_method = "rank1",
  control = list(maxiter = 1000),
  cc = 1e-3
)

normalized_structure_plot(log1p_mod_r1_init, loadings_order = 1:1000)

set.seed(1)
log1p_mod_r1_init_k5 <- fit_poisson_log1p_nmf(
  Y = Y,
  K = 5,
  loglik = "exact",
  init_method = "rank1",
  control = list(maxiter = 1000),
  cc = 1e-3
)

normalized_structure_plot(log1p_mod_r1_init_k5, loadings_order = 1:1000)


library(fastTopics)
set.seed(1)
ft_fit <- fit_poisson_nmf(
  X = Y, k = 4, init.method = "random", control = list(nc = 7),
  verbose = "none"
)

structure_plot(ft_fit, loadings_order = 1:1000)


set.seed(1)
log1p_mod_big <- fit_poisson_log1p_nmf(
  Y = Y[,951:1000],
  K = 4,
  loglik = "exact",
  init_method = "rank1",
  control = list(maxiter = 1000),
  cc = 1e-3
)

normalized_structure_plot(log1p_mod_big, loadings_order = 1:1000)


