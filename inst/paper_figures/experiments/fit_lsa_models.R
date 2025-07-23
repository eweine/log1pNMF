library(Matrix)
library(log1pNMF)

load("../data/panc_cyto_lsa.Rdata")

cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)

K <- 13

for (cc in cc_vec) {
  
  set.seed(1)
  fit <- fit_poisson_log1p_nmf(
    Y = counts,
    K = K,
    cc = cc,
    init_method = "rank1",
    loglik = "exact",
    control = list(maxiter = 250)
  )
  
  readr::write_rds(fit, glue::glue("log1p_k{K}_c_{cc}_lsa.rds"))
  
}

# finally, fit a topic model with a rank 1 initialization here
r1_fit <- fastTopics:::fit_pnmf_rank1(counts)

init_LL <- cbind(
  r1_fit$L,
  matrix(
    data = 1e-5,
    nrow = nrow(counts),
    ncol = K - 1
  )
)
rownames(init_LL) <- rownames(counts)

init_FF <- cbind(
  r1_fit$F,
  matrix(
    data = 1e-5,
    nrow = ncol(counts),
    ncol = K - 1
  )
)
rownames(init_FF) <- colnames(counts)

fit0 <- init_poisson_nmf(X = counts, L = init_LL, F = init_FF)

nmf_fit <- fit_poisson_nmf(
  X = counts,
  fit0 = fit0,
  numiter = 1000,
  control = list(nc = 7)
)
