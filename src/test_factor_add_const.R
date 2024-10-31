# here, I want to figure out if my new algorithm that fits alpha will
# lead to the log-likelihood being closer to that of the glmpca model.
# I expect that it will be, but I'm not sure if it may even exceed
# that log-likelihood because it has an additional parameter.

set.seed(10)
n <- 500
p <- 250
K <- 4

library(distr)
library(Matrix)
#library(MatrixExtra)

l_dist <- UnivarMixingDistribution(
  Dirac(0),
  Exp(rate = 3),
  mixCoeff = rep(1/2,2)
)

f_dist <- UnivarMixingDistribution(
  Dirac(0),
  Exp(rate = 3),
  mixCoeff = rep(1/2,2)
)

l_sampler <- distr::r(l_dist)
f_sampler <- distr::r(f_dist)

LL <- matrix(
  data = l_sampler(n = n * K),
  nrow = n,
  ncol = K
)

FF <- matrix(
  data = f_sampler(n = p * K),
  nrow = p,
  ncol = K
)

Lambda <- exp(LL %*% t(FF))

Y <- matrix(
  data = rpois(n = prod(dim(Lambda)), lambda = as.vector(Lambda)),
  nrow = n,
  ncol = p
)

Y_dense <- Y
Y <- as(Y, "CsparseMatrix")

ll_const <- sum(lfactorial(Y_dense))

cc <- 0.0000001

set.seed(1)
log1p_alpha <- passPCA:::fit_factor_model_log1p_exact_add_const(
  Y = Y,
  K = 4,
  maxiter = 100,
  cc = cc
)

H_alpha <- exp(log1p_alpha$U %*% t(log1p_alpha$V)) - cc

ll <- sum(
  dpois(
    x = as.vector(Y_dense),
    lambda = as.vector(H_alpha),
    log = TRUE
  )
)

ll_man <- tail(log1p_alpha$loglik)[1] - ll_const

library(fastglmpca)
set.seed(1)
glmpca_init <- init_glmpca_pois(
  col_size_factor = FALSE,
  row_intercept = FALSE,
  Y = Y,
  K = 4
)

glmpca_fit <- fit_glmpca_pois(
  Y = Y,
  fit0 = glmpca_init
)

ll2 <- sum(
  dpois(
    x = as.vector(Y_dense),
    lambda = as.vector(exp(glmpca_fit$U %*% diag(glmpca_fit$d) %*% t(glmpca_fit$V))),
    log = TRUE
  )
)

set.seed(1)
log1p <- passPCA::fit_factor_model_log1p_exact(
  Y = Y,
  K = 4,
  maxiter = 100000,
  s = rep(cc, n)
)
