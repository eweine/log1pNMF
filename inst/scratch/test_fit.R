# here, I want to test the matrix factorization code

library(passPCA)
library(Matrix)
set.seed(1)
n <- 1000
p <- 500
K <- 5

# first, simulate some data
dat <- generate_data_simple(n, p, K)

dat$U <- dat$U[rowSums(dat$Y) > 0, ]
dat$V <- dat$V[colSums(dat$Y) > 0, ]
dat$Y <- dat$Y[rowSums(dat$Y) > 0, ]
dat$Y <- dat$Y[, colSums(dat$Y) > 0]

library(rbenchmark)

library(tictoc)

tic()
fit <- fit_factor_model_log1p(dat$Y, K = 5, maxiter = 10)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

# These seem to match quite closely
#plot(as.vector(actual_lambda), as.vector(fitted_lambda))

tic()
fit_quad_approx <- fit_factor_model_log1p_quad_approx_full(
  dat$Y,
  K = 5,
  maxiter = 10,
  approx_range = c(0, 2)
)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

tic()
fit_lin_approx <- fit_factor_model_log1p_lin_approx_sparse(
  dat$Y,
  K = 5,
  maxiter = 10,
  a = 1
)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

