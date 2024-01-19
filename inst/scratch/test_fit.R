# here, I want to test the matrix factorization code

library(passPCA)
library(Matrix)
set.seed(1)
n <- 2500
p <- 1500
K <- 10

# first, simulate some data
dat <- generate_data_simple(n, p, K)

dat$U <- dat$U[rowSums(dat$Y) > 0, ]
dat$V <- dat$V[colSums(dat$Y) > 0, ]
dat$Y <- dat$Y[rowSums(dat$Y) > 0, ]
dat$Y <- dat$Y[, colSums(dat$Y) > 0]

fit <- fit_factor_model_log1p(dat$Y, K = 10, maxiter = 10)

actual_lambda <- exp(tcrossprod(dat$U, dat$V)) - 1
fitted_lambda <- exp(tcrossprod(fit$U, fit$V)) - 1

# These seem to match quite closely
plot(as.vector(actual_lambda), as.vector(fitted_lambda))

