# I think that the next thing I should try is to do more direct
# multiplication instead of using the link function
n_cells <- 400

grouping <- c(
  rep("A1", 100), rep("A2", 100),
  rep("B1", 100), rep("B2", 100)
)

set.seed(1)
lambda_group1_c <- c(rexp(200), rep(0, 800))
lambda_group2_c <- c(rep(0, 200), rexp(200), rep(0, 600))

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

lambda_list <- list()
Y_list <- list()

lambda_list[[1]] <- lambda_group1_c
lambda_list[[2]] <- lambda_group1_t
lambda_list[[3]] <- lambda_group2_c
lambda_list[[4]] <- lambda_group2_t

set.seed(1)
n_cells_per_group <- 100
n_genes <- 1000

for (group in 1:4) {
  
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
  matrix(data = 1e-5, nrow = nrow(Y), ncol = 2)
)

init_FF <- cbind(
  ft_r1$F,
  matrix(data = 1e-5, nrow = ncol(Y), ncol = 2)
)

ft_init <- init_poisson_nmf(X = Y, F = init_FF, L = init_LL)

ft_r1_init <- fit_poisson_nmf(
  X = Y, 
  fit0 = ft_init, 
  control = list(nc = 7),
  verbose = "none"
)
sp1 <- structure_plot(ft_r1_init, loadings_order = 1:n_cells, grouping = grouping, gap = 10) +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12)) + ylab("Membership") + ggtitle("Topic Model")

log1p_fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)

for (cc in cc_vec) {
  
  set.seed(1)
  log1p_fit_list[[as.character(cc)]] <- fit_poisson_log1p_nmf(
    Y = Y, K = 6, loglik = "exact", init_method = "rank1",
    control = list(maxiter = 250, verbose = FALSE), cc = cc
  )
  
}


sp2 <- normalized_structure_plot(log1p_fit_list[[as.character(1)]], loadings_order = 1:n_cells, grouping = grouping, gap = 10, topics = c(1, rev(2:6))) +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12)) + ylab("Membership") + ggtitle("log1p Model (c = 1)")

sp3 <- normalized_structure_plot(log1p_fit_list[[as.character(0.001)]], loadings_order = 1:n_cells, grouping = grouping, gap = 10, topics = c(1, rev(2:6))) +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12)) + ylab("Membership") + ggtitle("log1p Model (c = 1e-3)")


log1p_fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)

for (cc in cc_vec) {
  
  set.seed(1)
  log1p_fit_list[[as.character(cc)]] <- fit_poisson_log1p_nmf(
    Y = Y, K = 6, loglik = "exact", init_method = "rank1",
    control = list(maxiter = 250, verbose = FALSE), cc = cc
  )
  
}


sp2 <- normalized_structure_plot(log1p_fit_list[[as.character(1)]], loadings_order = 1:n_cells, grouping = grouping, gap = 10, topics = c(1, rev(2:6))) +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12)) + ylab("Membership") + ggtitle("log1p Model (c = 1)")

sp3 <- normalized_structure_plot(log1p_fit_list[[as.character(0.001)]], loadings_order = 1:n_cells, grouping = grouping, gap = 10, topics = c(1, rev(2:6))) +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12)) + ylab("Membership") + ggtitle("log1p Model (c = 1e-3)")
