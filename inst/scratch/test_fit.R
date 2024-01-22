# here, I want to test the matrix factorization code

library(passPCA)
library(Matrix)
set.seed(1)
n <- 1500
p <- 750
K <- 6

# first, simulate some data
dat <- generate_data_simple(n, p, K)

dat$U <- dat$U[rowSums(dat$Y) > 0, ]
dat$V <- dat$V[colSums(dat$Y) > 0, ]
dat$Y <- dat$Y[rowSums(dat$Y) > 0, ]
dat$Y <- dat$Y[, colSums(dat$Y) > 0]

library(rbenchmark)

library(tictoc)

tic()
fit <- fit_factor_model_log1p(dat$Y, K = 6, maxiter = 10)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

print(cor(as.vector(actual_lambda), as.vector(fitted_lambda)))

# These seem to match quite closely
#plot(as.vector(actual_lambda), as.vector(fitted_lambda))

tic()
fit <- fit_factor_model_log1p_quad_approx_full(
  dat$Y,
  K = 6,
  maxiter = 10,
  approx_range = c(0, 2)
)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

print(cor(as.vector(actual_lambda), as.vector(fitted_lambda)))

tic()
fit <- fit_factor_model_log1p_quad_approx_sparse(
  dat$Y,
  K = 6,
  maxiter = 10,
  approx_range = c(0, 2)
)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

print(cor(as.vector(actual_lambda), as.vector(fitted_lambda)))


tic()
fit <- fit_factor_model_log1p_lin_approx_sparse(
  dat$Y,
  K = 6,
  maxiter = 10,
  a = 1
)
toc()

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

print(cor(as.vector(actual_lambda), as.vector(fitted_lambda)))
