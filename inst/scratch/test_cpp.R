# generate data
set.seed(1)
n <- 500
p <- 5

X <- matrix(
  data = 0, nrow = n, ncol = p
)

X[, 1] <- abs(rnorm(n, sd = .5))
X[, 2] <- abs(rnorm(n, sd = 1))
X[, 3] <- rexp(n, rate = 3)
X[, 4] <- rgamma(n, shape = 3, rate = 3)
X[1:460, 5] <- abs(rnorm(460, sd = .01))
X[461:500, 5] <- rexp(40, 1)

b <- c(0, 0, .1, .3, .6)

lambda <- exp(X %*% b) - 1
y <- rpois(n, lambda)

set.seed(1)
init <- abs(rnorm(5))

# first, I want to optimize there in R

get_sparse_term_loglik <- function(y_nz, X_nz, b) {

  sum(
    y_nz * log(exp(X_nz %*% b) - 1)
  )

}

get_exp_term_loglik_exact <- function(X, b) {

  -sum(exp(X %*% b))

}

get_exp_term_loglik_lin_approx <- function(X, b, a1) {

  -a1 * sum(colSums(X) * b)

}

get_exp_term_loglik_quad_approx <- function(X, b, a1, a2) {

  get_exp_term_loglik_lin_approx(X, b, a1) - a2 * (
    t(b) %*% crossprod(X) %*% b
  )

}

# exact log-likelihood for regression
get_exact_loglik <- function(b, X, y) {

  y_nz_idx <- which(y != 0)
  ll <- get_sparse_term_loglik(
    y[y_nz_idx], X[y_nz_idx, , drop = FALSE], b
  ) + get_exp_term_loglik_exact(X, b)
  return(ll)

}

# log-likelihood approximating only sparse terms with linear approximation
get_lin_sparse_approx_loglik <- function(b, X, y, a1) {

  y_nz_idx <- which(y != 0)
  y_z_idx <- which(y == 0)
  ll <- get_sparse_term_loglik(
    y[y_nz_idx], X[y_nz_idx, , drop = FALSE], b
  ) + get_exp_term_loglik_exact(
    X[y_nz_idx, , drop = FALSE], b
  ) + get_exp_term_loglik_lin_approx(
    X[y_z_idx, , drop = FALSE], b, a1
  )
  return(ll)

}

# log-likelihood approximating only sparse terms with quadratic approximation
get_quad_sparse_approx_loglik <- function(b, X, y, a1, a2) {

  y_nz_idx <- which(y != 0)
  y_z_idx <- which(y == 0)
  ll <- get_sparse_term_loglik(
    y[y_nz_idx], X[y_nz_idx, , drop = FALSE], b
  ) + get_exp_term_loglik_exact(
    X[y_nz_idx, , drop = FALSE], b
  ) + get_exp_term_loglik_quad_approx(
    X[y_z_idx, , drop = FALSE], b, a1, a2
  )
  return(ll)

}

# log-likelihood approximating ALL terms with quadratic approximation
get_quad_full_approx_loglik <- function(b, X, y, a1, a2) {

  y_nz_idx <- which(y != 0)
  ll <- get_sparse_term_loglik(
    y[y_nz_idx], X[y_nz_idx, , drop = FALSE], b
  ) + get_exp_term_loglik_quad_approx(
    X, b, a1, a2
  )
  return(ll)

}

exact_sol <- optim(
  par = init,
  fn = get_exact_loglik,
  X = X,
  y = y,
  lower = rep(0, 5),
  control = list(fnscale = -1),
  method = "L-BFGS-B"
)$par

# now, need to determine the lengths of the approximation intervals
approx_poly <- pracma::polyApprox(
  exp,
  a = 0,
  b = log(5),
  n = 1
)

a0_linear <- approx_poly$p[2]
a1_linear <- approx_poly$p[1]

lin_sparse_approx_sol <- optim(
  par = init,
  fn = get_lin_sparse_approx_loglik,
  X = X,
  y = y,
  a1 = a1_linear,
  lower = rep(0, 5),
  control = list(fnscale = -1),
  method = "L-BFGS-B"
)$par


approx_poly <- pracma::polyApprox(
  exp,
  a = 0,
  b = log(5),
  n = 2
)

a0_quad <- approx_poly$p[3]
a1_quad <- approx_poly$p[2]
a2_quad <- approx_poly$p[1]

quad_sparse_approx_sol <- optim(
  par = init,
  fn = get_quad_sparse_approx_loglik,
  X = X,
  y = y,
  a1 = a1_quad,
  a2 = a2_quad,
  lower = rep(0, 5),
  control = list(fnscale = -1),
  method = "L-BFGS-B"
)$par

approx_poly <- pracma::polyApprox(
  exp,
  a = 0,
  b = log(8),
  n = 2
)

a1_quad_full <- approx_poly$p[2]
a2_quad_full <- approx_poly$p[1]

quad_full_approx_sol <- optim(
  par = init,
  fn = get_quad_full_approx_loglik,
  X = X,
  y = y,
  a1 = a1_quad_full,
  a2 = a2_quad_full,
  lower = rep(0, 5),
  control = list(fnscale = -1),
  method = "L-BFGS-B"
)$par

y_nz_idx <- which(y != 0)
y_z_idx <- which(y == 0)

library(tictoc)

set.seed(1)
init <- abs(rnorm(5))
# for the exact method

tic()
b_cpp_exact <- passPCA:::solve_pois_reg_log1p (
  X,
  y[y != 0],
  as.integer(which(y != 0) - 1),
  init,
  0:4,
  200,
  .01,
  .25
)
toc()

set.seed(1)
init <- abs(rnorm(5))

tic()
b_cpp_linear <-  passPCA:::solve_pois_reg_log1p_lin_approx (
  X[y_nz_idx, ],
  y[y != 0],
  a1_linear * colSums(X[y_z_idx, ]),
  init,
  0:4,
  200,
  .01,
  .25
)
toc()

set.seed(1)
init <- abs(rnorm(5))

tic()
b_cpp_quad <-  passPCA:::solve_pois_reg_log1p_quad_approx (
  X[y_nz_idx, ],
  y[y != 0],
  a1_quad * colSums(X[y_z_idx, ]),
  crossprod(X[y_z_idx, ]),
  a2_quad,
  init,
  0:4,
  200,
  .01,
  .25
)
toc()

set.seed(1)
init <- abs(rnorm(5))

tic()
b_cpp_quad_full <-  passPCA:::solve_pois_reg_log1p_quad_approx_full (
  X[y_nz_idx, ],
  y[y != 0],
  a1_quad_full * colSums(X),
  crossprod(X),
  a2_quad_full,
  init,
  0:4,
  200,
  .01,
  .25
)
toc()

# everything seems to be agreeing nearly exactly

# the last glm method I need to implement is the exact method.
# once I have done this, then I will have everything needed for the
# factor model to work.

# Now, I want to implement a factor model pipeline.
# I will start with the exact method,
# then move to the full quadratic approximation,
# and then I will try and figure out the other two methods
# as they will require a bit of work

