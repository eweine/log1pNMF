library(Matrix)
set.seed(5)
n_genes <- 500
n_cells <- 2500
K <- 4
library(passPCA)
FF <- matrix(
  data = 0,
  nrow = n_genes,
  ncol = K
)

LL <- matrix(
  data = 0,
  nrow = n_cells,
  ncol = K
)

FF[, 1] <- c(rep(log(1.00225), 250), rep(log(1.5), 250))
FF[210:260, 2] <- 1.8
FF[250:300, 3] <- 1.8
FF[290:340, 4] <- 1.8

LL[, 1] <- 1

for (j in 2:K) {

  LL[, j] <- sample(c(0, 1), size = n_cells, replace = TRUE, prob = c(0.8, 0.2))

}

H <- tcrossprod(LL, FF)
Lambda <- exp(H) - 1

y <- rpois(n_cells * n_genes, lambda = as.vector(Lambda))
Y <- matrix(
  data = y,
  nrow = n_cells,
  ncol = n_genes
)
Y <- as(Y, "sparseMatrix")

library(fastTopics)

ft_mod <- fit_poisson_nmf(Y, k=4, numiter = 250)

log1p_mod <- fit_factor_model_log1p_quad_approx_sparse(
  Y=Y, K = 4, approx_range = c(0, 1.25), maxiter = 250
)

true_factor_df <- data.frame(
  val = c(
    FF[,1], FF[,2], FF[,3], FF[,4]
  ),
  factor = c(
    rep("1", 500), rep("2", 500), rep("3", 500), rep("4", 500)
  ),
  pos = c(rep(1:500, 4))
)
library(ggplot2)
ggplot(true_factor_df, aes(x = pos, y = val)) +
  geom_bar(stat = "identity") +  # Use bars to represent lambda values
  facet_wrap(~ factor) +          # Create a panel for each group
  labs(x = "Position", y = "Value") +
  cowplot::theme_cowplot()  # Use a minimal theme for a clean look

ft_factor_df <- data.frame(
  val = c(
    ft_mod$F[,1], ft_mod$F[,2], ft_mod$F[,3], ft_mod$F[,4]
  ),
  factor = c(
    rep("1", 500), rep("2", 500), rep("3", 500), rep("4", 500)
  ),
  pos = c(rep(1:500, 4))
)
library(ggplot2)
ggplot(ft_factor_df, aes(x = pos, y = val)) +
  geom_bar(stat = "identity") +  # Use bars to represent lambda values
  facet_wrap(~ factor) +          # Create a panel for each group
  labs(x = "Position", y = "Value") +
  cowplot::theme_cowplot()  # Use a minimal theme for a clean look


max_col <- apply(log1p_mod$U, 2, max)
log1p_LL <- sweep(log1p_mod$U, 2, max_col, FUN = "/")
log1p_FF <- sweep(log1p_mod$V, 2, max_col, FUN = "*")

log1p_factor_df <- data.frame(
  val = c(
    log1p_FF[,1], log1p_FF[,2], log1p_FF[,3], log1p_FF[,4]
  ),
  factor = c(
    rep("1", 500), rep("2", 500), rep("3", 500), rep("4", 500)
  ),
  pos = c(rep(1:500, 4))
)
library(ggplot2)
ggplot(log1p_factor_df, aes(x = pos, y = val)) +
  geom_bar(stat = "identity") +  # Use bars to represent lambda values
  facet_wrap(~ factor) +          # Create a panel for each group
  labs(x = "Position", y = "Value") +
  cowplot::theme_cowplot()  # Use a minimal theme for a clean look


# I think this is a pretty good demonstration of what
# I am intending to show.
# I think in addition to structure plots, this provides a
# good demonstration.
