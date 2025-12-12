library(Matrix)
library(log1pNMF)
library(dplyr)
library(NNLM)

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")

set.seed(1)

s <- Matrix::rowSums(counts) 
s <- s / mean(s)
Y_tilde <- Matrix::Diagonal(x = 1/s) %*% counts
Y_tilde <- MatrixExtra::mapSparse(Y_tilde, log1p)

set.seed(1)
fit0 <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 13,
  init_method = "rank1",
  loglik = "exact",
  control = list(maxiter = 1)
)

fit_approx <- fit_poisson_log1p_nmf(
  Y = counts,
  K = 13,
  loglik = "approx",
  control = list(maxiter = 249),
  init_LL = fit0$LL,
  init_FF = fit0$FF
)
readr::write_rds(fit_approx, "~/Documents/data/passPCA/lsa_k13_cheby_approx.rds")

init_W <- fit0$LL
init_H <- t(fit0$FF)

fit <- nnmf(
  A = as.matrix(Y_tilde),
  init = list(W = init_W, H = init_H),
  k = 13
)

readr::write_rds(fit, "~/Documents/data/passPCA/lsa_k13_frob_fit.rds")
