# I think that the next thing I should try is to do more direct
# multiplication instead of using the link function

set.seed(1)
lambda_group1_c <- c(rexp(200), rep(0, 800))
lambda_group2_c <- c(rep(0, 200), rexp(200), rep(0, 600))
lambda_group3_c <- c(rep(0, 400), rexp(200), rep(0, 400))
lambda_group4_c <- c(rep(0, 600), rexp(200), rep(0, 200))
lambda_group5_c <- c(rep(0, 800), rexp(200))

mult_vec <- rep(2.5, 1000)
mult_vec[sample(1:1000, size = 800)] <- 0

lambda_group1_t <- if_else(
  mult_vec == 0,
  lambda_group1_c,
  exp(mult_vec + lambda_group1_c)
)

lambda_group2_t <- if_else(
  mult_vec == 0,
  lambda_group2_c,
  exp(mult_vec + lambda_group2_c)
)

lambda_group3_t <- if_else(
  mult_vec == 0,
  lambda_group3_c,
  exp(mult_vec + lambda_group3_c)
)

lambda_group4_t <- if_else(
  mult_vec == 0,
  lambda_group4_c,
  exp(mult_vec + lambda_group4_c)
)

lambda_group5_t <- if_else(
  mult_vec == 0,
  lambda_group4_c,
  exp(mult_vec + lambda_group5_c)
)

lambda_list <- list()
Y_list <- list()

lambda_list[[1]] <- lambda_group1_c
lambda_list[[2]] <- lambda_group1_t
lambda_list[[3]] <- lambda_group2_c
lambda_list[[4]] <- lambda_group2_t
lambda_list[[5]] <- lambda_group3_c
lambda_list[[6]] <- lambda_group3_t
lambda_list[[7]] <- lambda_group4_c
lambda_list[[8]] <- lambda_group4_t
lambda_list[[9]] <- lambda_group5_c
lambda_list[[10]] <- lambda_group5_t

set.seed(1)
n_cells_per_group <- 100
n_genes <- 1000

for (group in 1:10) {
  
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

Y <- as(Y, "CsparseMatrix")
Y <- Y[,Matrix::colSums(Y) > 0]

library(fastTopics)

ft_r1 <- fastTopics:::fit_pnmf_rank1(Y)

init_LL <- cbind(
  ft_r1$L,
  matrix(data = 1e-5, nrow = nrow(Y), ncol = 5)
)

init_FF <- cbind(
  ft_r1$F,
  matrix(data = 1e-5, nrow = ncol(Y), ncol = 5)
)

ft_init <- init_poisson_nmf(X = Y, F = init_FF, L = init_LL)

ft_r1_init <- fit_poisson_nmf(X = Y, fit0 = ft_init)
structure_plot(ft_r1_init, loadings_order = 1:n_cells)



set.seed(1)
log1p_fit <- fit_poisson_log1p_nmf(
  Y = Y, K = 6, loglik = "exact", init_method = "rank1",
  control = list(maxiter = 250), cc = 1
)

normalized_structure_plot(log1p_fit, loadings_order = 1:n_cells)


