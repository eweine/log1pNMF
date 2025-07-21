# here, I want to run a hierarchical simulation to understand the difference
# between the topic model and log1p

library(fastTopics)
library(log1pNMF)
set.seed(1)

LL <- matrix(
  data = c(
    rep(1, 300), rep(0, 600), # k1
    rep(0, 300), rep(1, 300), rep(0, 300), # k2
    rep(0, 600), rep(1, 300), # k3
    rep(1, 200), rep(0, 700), # k4
    rep(0, 300), rep(1, 200), rep(0, 400), # k5
    rep(0, 600), rep(1, 200), rep(0, 100), # k6
    rep(1, 100), rep(0, 800), # k7
    rep(0, 300), rep(1, 100), rep(0, 500), # k8
    rep(0, 600), rep(1, 100), rep(0, 200) # k9
  ),
  nrow = 900, 
  ncol = 9
)

# if FF is orthogonal then additive vs. multiplicative shouldn't actually matter
# which could either be viewed as a feature or a bug
FF <- matrix(
  data = c(
    pmax(rexp(100), 5e-2), rep(0, 800),
    rep(0, 100), pmax(rexp(100), 5e-2), rep(0, 700),
    rep(0, 200), pmax(rexp(100), 5e-2), rep(0, 600),
    rep(0, 300), pmax(rexp(100), 5e-2), rep(0, 500),
    rep(0, 400), pmax(rexp(100), 5e-2), rep(0, 400),
    rep(0, 500), pmax(rexp(100), 5e-2), rep(0, 300),
    rep(0, 600), pmax(rexp(100), 5e-2), rep(0, 200),
    rep(0, 700), pmax(rexp(100), 5e-2), rep(0, 100),
    rep(0, 800), pmax(rexp(100), 5e-2)
  ),
  nrow = 900,
  ncol = 9
)

B <- tcrossprod(LL, FF)

Lambda <- expm1(B)

Y <- matrix(
  data = rpois(900 * 900, lambda = as.vector(Lambda)),
  nrow = 900,
  ncol = 900
)

Y <- as(Y, "CsparseMatrix")
set.seed(1)
ft_fit <- fit_poisson_nmf(X = Y, k = 9, numiter = 1000, control = list(nc = 7))

structure_plot(LL, loadings_order = 1:900)
structure_plot(ft_fit$L, loadings_order = 1:900)

set.seed(1)
log1p_fit <- fit_poisson_log1p_nmf(
  Y = Y, K = 9, loglik = "exact", init_method = "random",
  control = list(maxiter = 1000)
)

structure_plot(log1p_fit$LL, loadings_order = 1:900)

#set.seed(1)
#log1p_fit_r1_init <- fit_poisson_log1p_nmf(
#  Y = Y, K = 9, loglik = "exact", init_method = "rank1",
#  control = list(maxiter = 1000)
#)

#normalized_structure_plot(log1p_fit_r1_init, loadings_order = 1:900)
