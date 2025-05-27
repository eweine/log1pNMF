library(distr)

set.seed(10)
n <- 10000
p <- 3
cc <- 3.0

X_dist <- UnivarMixingDistribution(
  Exp(rate = 3),
  Unif(1, 1.025),
  mixCoeff = c((3/5), (2/5))
)

X_sampler <- distr::r(X_dist)

X <- matrix(
  data = X_sampler(n * p),
  nrow = n,
  ncol = p
)

b <- c(0.1, 0.25, .5)
s <- runif(n, 0.45, 1.55)

lambda <- s * cc * (exp(X %*% b) - 1)

y <- rpois(n = n, lambda = lambda)

# now, I would like to build code that can estimate both b and c
# I want to start with using CCD and update c at each iteration
y_nz_R_idx <- which(y != 0)

poly_approx <- pracma::polyApprox(
  exp,
  0,
  1,
  2
)

a0 <- poly_approx$p[3]
a1 <- poly_approx$p[2]
a2 <- poly_approx$p[1]

fit_out <- log1pNMF:::solve_pois_reg_log1p_quad_approx_sparse_vec_s_fit_c(
  t(X),
  y[y_nz_R_idx],
  y_nz_R_idx - 1,
  sum(y),
  s[y_nz_R_idx],
  sum(s),
  sum(s) - sum(s[y_nz_R_idx]),
  crossprod(X, s),
  t(X) %*% diag(s) %*% X,
  a0,
  a1,
  a2,
  c(0, .01, 0.1),
  3.0,
  0:2,
  200,
  .01,
  .25
)

# Now, it seems like I'm able to successfully fit b when there is a value of c
# the next question is if I am able to fit c well

lambda_fit <- s * fit_out$c * (exp(X %*% fit_out$b) - 1)

plot(y, as.vector(lambda_fit))
