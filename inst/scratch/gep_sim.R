# I think to start I should attempt to see how the model might perform 
# when there is overdispersion.

n_genes <- 1000
n_groups <- 4
n_cells_per_group <- 250

high_expr <- 25
low_expr <- 1

lambda_list <- list()
Y_list <- list()

lambda_list[[1]] <- c(
  rep(high_expr, 250),
  rep(0, 500),
  rep(low_expr, 250)
)

lambda_list[[2]] <- c(
  rep(0, 250),
  rep(high_expr, 250),
  rep(low_expr, 250),
  rep(0, 250)
)

lambda_list[[3]] <- c(
  rep(0, 250),
  rep(low_expr, 250),
  rep(high_expr, 250),
  rep(0, 250)
)

lambda_list[[4]] <- c(
  rep(low_expr, 250),
  rep(0, 250),
  rep(0, 250),
  rep(high_expr, 250)
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

library(fastTopics)

set.seed(1)
ft_fit <- fit_poisson_nmf(X = Y, k = 4, init.method = "random")

structure_plot(ft_fit, loadings_order = 1:1000)

library(log1pNMF)

set.seed(1)
log1p_mod <- fit_poisson_log1p_nmf(
  Y = Y,
  K = 4,
  loglik = "exact",
  init_method = "rank1",
  control = list(maxiter = 250),
  cc = 1e-3
)

normalized_structure_plot(log1p_mod, loadings_order = 1:1000)