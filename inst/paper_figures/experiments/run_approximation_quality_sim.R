n <- 500
p <- 500
K <- 5
library(Matrix)
library(log1pNMF)
library(NNLM)

# I would like the simulation to result in fairly sparse data...
library(distr)
LL <- matrix(nrow = n,ncol = K)
FF <- matrix(nrow = p,ncol = K)
set.seed(1)

for (k in 1:K) {
  
  rate_l <- runif(1, 0.15, 5)
  l_dist <- UnivarMixingDistribution(
    Unif(1e-12,0.1),
    Exp(rate_l),
    mixCoeff = c(0.9, 0.1)
  )
  
  rate_f <- runif(1, 0.15, 5)
  f_dist <- UnivarMixingDistribution(
    Unif(1e-12,0.1),
    Exp(rate_f),
    mixCoeff = c(0.9, 0.1)
  )
  
  l_sampler <- distr::r(l_dist)
  f_sampler <- distr::r(f_dist)
  
  LL[,k] <- l_sampler(n)
  FF[,k] <- f_sampler(p)
  
}

Lambda <- tcrossprod(LL, FF)

Y <- matrix(
  data = rpois(n * p, as.vector(Lambda)),
  nrow = n, ncol = p
)

Y <- as(Y, "CsparseMatrix")
cc_vec <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4)

fit_list_exact <- list()
fit_list_approx <- list()

for (cc in cc_vec) {
  
  set.seed(1)
  fit_exact <- fit_poisson_log1p_nmf(
    Y = Y,
    K = K,
    cc = cc,
    s = FALSE,
    loglik = "exact",
    init_method = "random",
    control = list(
      maxiter = 1e5,
      threads = 7,
      tol = 1e-6
    )
  )
  
  set.seed(1)
  fit_approx <- fit_poisson_log1p_nmf(
    Y = Y,
    K = K,
    cc = cc,
    s = FALSE,
    loglik = "approx",
    init_method = "random",
    control = list(
      maxiter = 1e5,
      threads = 7,
      tol = 1e-6
    )
  )
  
  Y_tilde <- MatrixExtra::mapSparse(Y, function(x){log1p(x/cc)})
  frob_fit <- NNLM::nnmf(A = as.matrix(Y_tilde), k = K)
  
  fit_list_exact[[as.character(cc)]] <- fit_exact
  fit_list_approx[[as.character(cc)]] <- fit_approx
  
}

relative_rmse <- function(A, B) {
  # Ensure matrices have the same dimensions
  if (!all(dim(A) == dim(B))) {
    stop("Matrices A and B must have the same dimensions.")
  }
  
  # Compute relative RMSE
  sqrt(mean((A - B)^2)) / sqrt(mean(B^2))
}

for (cc in cc_vec) {
  
  fit_list_exact[[as.character(cc)]]$ll <- sum(
    dpois(
      as.vector(as.matrix(Y)),
      as.vector(fitted(fit_list_exact[[as.character(cc)]])),
      log = TRUE
    )
  ) / (n * p)
  fit_list_approx[[as.character(cc)]]$ll <- sum(
    dpois(
      as.vector(as.matrix(Y)),
      as.vector(fitted(fit_list_approx[[as.character(cc)]])),
      log = TRUE
    )
  ) / (n * p)
  fit_list_frob[[as.character(cc)]]$ll <- sum(
    dpois(
      as.vector(as.matrix(Y)),
      cc * (exp(fit_list_frob[[as.character(cc)]]$W %*% fit_list_frob[[as.character(cc)]]$H) - 1),
      log = TRUE
    )
  ) / (n * p)
  
  fit_list_exact[[as.character(cc)]]$rrmse <- relative_rmse(fitted(fit_list_exact[[as.character(cc)]]), Lambda)
  fit_list_approx[[as.character(cc)]]$rrmse <- relative_rmse(fitted(fit_list_approx[[as.character(cc)]]), Lambda)
  fit_list_frob[[as.character(cc)]]$rrmse <- relative_rmse(cc * (exp(fit_list_frob[[as.character(cc)]]$W %*% fit_list_frob[[as.character(cc)]]$H) - 1), Lambda)
  
}

save(fit_list_exact, fit_list_approx, file = "~/Documents/data/fit_list_sim_tm.Rdata")

n <- 500
p <- 500
K <- 5
library(Matrix)
library(log1pNMF)
# I would like the simulation to result in fairly sparse data...
library(distr)
LL <- matrix(nrow = n,ncol = K)
FF <- matrix(nrow = p,ncol = K)
set.seed(1)

for (k in 1:K) {
  
  rate_l <- runif(1, 1.25, 10)
  l_dist <- UnivarMixingDistribution(
    Unif(1e-12,0.16),
    Exp(rate_l),
    mixCoeff = c(0.9, 0.1)
  )
  
  rate_f <- runif(1, 1.25, 10)
  f_dist <- UnivarMixingDistribution(
    Unif(1e-12,0.16),
    Exp(rate_f),
    mixCoeff = c(0.9, 0.1)
  )
  
  l_sampler <- distr::r(l_dist)
  f_sampler <- distr::r(f_dist)
  
  LL[,k] <- l_sampler(n)
  FF[,k] <- f_sampler(p)
  
}

Lambda <- exp(tcrossprod(LL, FF)) - 1

Y <- matrix(
  data = rpois(n * p, as.vector(Lambda)),
  nrow = n, ncol = p
)

Y <- as(Y, "CsparseMatrix")
length(Y@x) / prod(dim(Y))

cc_vec <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4)

fit_list_exact <- list()
fit_list_approx <- list()

for (cc in cc_vec) {
  
  set.seed(1)
  fit_exact <- fit_poisson_log1p_nmf(
    Y = Y,
    K = K,
    cc = cc,
    s = FALSE,
    loglik = "exact",
    init_method = "random",
    control = list(
      maxiter = 1e5,
      threads = 7,
      tol = 1e-6
    )
  )
  
  set.seed(1)
  fit_approx <- fit_poisson_log1p_nmf(
    Y = Y,
    K = K,
    cc = cc,
    s = FALSE,
    loglik = "approx",
    init_method = "random",
    control = list(
      maxiter = 1e5,
      threads = 7,
      tol = 1e-6
    )
  )
  
  fit_list_exact[[as.character(cc)]] <- fit_exact
  fit_list_approx[[as.character(cc)]] <- fit_approx
  
}


for (cc in cc_vec) {
  
  fit_list_exact[[as.character(cc)]]$ll <- sum(
    dpois(
      as.vector(as.matrix(Y)),
      as.vector(fitted(fit_list_exact[[as.character(cc)]])),
      log = TRUE
    )
  ) / (n * p)
  fit_list_approx[[as.character(cc)]]$ll <- sum(
    dpois(
      as.vector(as.matrix(Y)),
      as.vector(fitted(fit_list_approx[[as.character(cc)]])),
      log = TRUE
    )
  ) / (n * p)
  
}

save(fit_list_exact, fit_list_approx, file = "~/Documents/data/fit_list_sim_c1.Rdata")

n <- 500
p <- 500
K <- 5
library(Matrix)
library(log1pNMF)

# I would like the simulation to result in fairly sparse data...
library(distr)
LL <- matrix(nrow = n,ncol = K)
FF <- matrix(nrow = p,ncol = K)
set.seed(1)

for (k in 1:K) {
  
  rate_l <- runif(1, 1.2, 5)
  l_dist <- UnivarMixingDistribution(
    Unif(0.35,1.43),
    Exp(rate_l),
    mixCoeff = c(0.9, 0.1)
  )
  
  rate_f <- runif(1, 1.2, 5)
  f_dist <- UnivarMixingDistribution(
    Unif(0.35,1.43),
    Exp(rate_f),
    mixCoeff = c(0.9, 0.1)
  )
  
  l_sampler <- distr::r(l_dist)
  f_sampler <- distr::r(f_dist)
  
  LL[,k] <- l_sampler(n)
  FF[,k] <- f_sampler(p)
  
}

Lambda <- 1e-3 * (exp(tcrossprod(LL, FF)) - 1)

Y <- matrix(
  data = rpois(n * p, as.vector(Lambda)),
  nrow = n, ncol = p
)

Y <- as(Y, "CsparseMatrix")
length(Y@x) / prod(dim(Y))

cc_vec <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4)

fit_list_exact <- list()
fit_list_approx <- list()

for (cc in cc_vec) {
  
  set.seed(1)
  fit_exact <- fit_poisson_log1p_nmf(
    Y = Y,
    K = K,
    cc = cc,
    s = FALSE,
    loglik = "exact",
    init_method = "random",
    control = list(
      maxiter = 1e5,
      threads = 7,
      tol = 1e-6
    )
  )
  
  set.seed(1)
  fit_approx <- fit_poisson_log1p_nmf(
    Y = Y,
    K = K,
    cc = cc,
    s = FALSE,
    loglik = "approx",
    init_method = "random",
    control = list(
      maxiter = 1e5,
      threads = 7,
      tol = 1e-6
    )
  )
  
  fit_list_exact[[as.character(cc)]] <- fit_exact
  fit_list_approx[[as.character(cc)]] <- fit_approx
  
}


for (cc in cc_vec) {
  
  fit_list_exact[[as.character(cc)]]$ll <- sum(
    dpois(
      as.vector(as.matrix(Y)),
      as.vector(fitted(fit_list_exact[[as.character(cc)]])),
      log = TRUE
    )
  ) / (n * p)
  fit_list_approx[[as.character(cc)]]$ll <- sum(
    dpois(
      as.vector(as.matrix(Y)),
      as.vector(fitted(fit_list_approx[[as.character(cc)]])),
      log = TRUE
    )
  ) / (n * p)
  
}

save(fit_list_exact, fit_list_approx, file = "~/Documents/data/fit_list_sim_c1e-3.Rdata")

