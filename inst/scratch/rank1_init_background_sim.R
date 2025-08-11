# I think to start I should attempt to see how the model might perform 
# when there is overdispersion.

n_genes <- 1000
n_groups <- 5
n_cells_per_group <- 200
high_genes_per_group <- 500

high_expr <- 15
low_expr <- 1.75

lambda_list <- list()
Y_list <- list()

set.seed(1)

for (group in 1:n_groups) {
  
  lambda <- rep(low_expr, n_genes)
  high_idx <- sample(1:n_genes, size = high_genes_per_group)
  lambda[high_idx] <- high_expr
  lambda_list[[group]] <- lambda
  Lambda <- matrix(
    data = rep(lambda, n_cells_per_group),
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
ft_fit <- fit_poisson_nmf(
  X = Y, 
  k = 5, 
  init.method = "random", 
  control = list(nc = 7)
  )

structure_plot(ft_fit, loadings_order = 1:1000)

library(log1pNMF)

set.seed(1)
log1p_mod <- fit_poisson_log1p_nmf(
  Y = Y,
  K = 5,
  loglik = "exact",
  init_method = "rank1",
  control = list(maxiter = 250),
  cc = 1
)

normalized_structure_plot(log1p_mod, loadings_order = 1:1000)
