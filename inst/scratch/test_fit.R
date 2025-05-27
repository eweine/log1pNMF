# here, I want to test the matrix factorization code

library(log1pNMF)
library(Matrix)
set.seed(1)
n <- 2500
p <- 1250
K <- 5

# so, now I've actually generated the data so that it is
# somewhat representative of real single cell data.
# Now, it would be good to look at the loglik vs. iter
# for these different models
# I have to edit the code to calculate the real log-likelihood

#d <- simulate_poisson_gene_data(n, p, K, 1, .5, sparse = TRUE)

#dat$Y <- dat$Y[rowSums(dat$Y) > 0, ]
#dat$Y <- dat$Y[, colSums(dat$Y) > 0]

# first, simulate some data
dat <- generate_data_simple(n, p, K)

dat$U <- dat$U[rowSums(dat$Y) > 0, ]
dat$V <- dat$V[colSums(dat$Y) > 0, ]
dat$Y <- dat$Y[rowSums(dat$Y) > 0, ]
dat$Y <- dat$Y[, colSums(dat$Y) > 0]

library(rbenchmark)

library(tictoc)

set.seed(1)
tic()
fit <- fit_factor_model_log1p(dat$Y, K = 5, maxiter = 50)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

# These seem to match quite closely
#plot(as.vector(actual_lambda), as.vector(fitted_lambda))
set.seed(1)

tic()
fit_quad_approx <- fit_factor_model_log1p_quad_approx_full(
  dat$Y,
  K = 5,
  maxiter = 50,
  approx_range = c(0, 2)
)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

set.seed(1)
tic()
fit_quad_approx_sparse <- fit_factor_model_log1p_quad_approx_sparse(
  dat$Y,
  K = 5,
  maxiter = 50,
  approx_range = c(0, .75)
)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

set.seed(1)
tic()
fit_lin_approx <- fit_factor_model_log1p_lin_approx_sparse(
  dat$Y,
  K = 5,
  maxiter = 50,
  a = 1
)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1


out_df <- data.frame(
  objective = c(
    fit$loglik,
    fit_quad_approx$loglik_exact,
    fit_quad_approx_sparse$loglik_exact,
    fit_lin_approx$loglik_exact
  ),
  iter = rep(0:50, 4),
  algorithm = c(
    rep("Exact", 51),
    rep("Full Quad Approx", 51),
    rep("Sparse Quad Approx", 51),
    rep("Sparse Lin Approx", 51)
  )
)

out_df$objective <- out_df$objective - sum(MatrixExtra::mapSparse(dat$Y, lfactorial))

library(dplyr)

out_df <- out_df %>%
  filter(algorithm != "Full Quad Approx")

library(ggplot2)

ggplot(data = out_df) +
  geom_line(aes(x = iter, y = objective, color = algorithm)) +
  cowplot::theme_cowplot() +
  ylab("Log-Likelihood") +
  xlab("Iteration")

out_df <- out_df %>%
  dplyr::mutate(
    diff_from_best = abs(objective - max(objective)) + 1
  )


ggplot(data = out_df) +
  geom_line(aes(x = iter, y = diff_from_best, color = algorithm)) +
  cowplot::theme_cowplot() +
  ylab("Log-Likelihood") +
  xlab("Iteration") +
  scale_y_continuous(trans = "log10")

# I think it would also be useful to plot fitted vs. mle lambdas
mle_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1
quad_approx_lambda <- exp(tcrossprod(fit_quad_approx$U, fit_quad_approx$V)) - 1
lin_approx_lambda <- exp(tcrossprod(fit_lin_approx$U, fit_lin_approx$V)) - 1

lambda_df <- data.frame(
  mle = as.vector(mle_lambda),
  quad_approx = as.vector(quad_approx_lambda),
  lin_approx = as.vector(lin_approx_lambda)
)

# To be honest, it would probably be easier to visualize
# the 0s and non-zero separately

ggplot(lambda_df) +
  geom_point(aes(x = mle, y = lin_approx), alpha = .3)



