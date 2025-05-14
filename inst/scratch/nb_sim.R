library(passPCA)
set.seed(1)
n <- 400
p <- 400
Y <- matrix(
  data = rnbinom(n = n * p, mu = 1, size = 0.1),
  nrow= n,
  ncol = p
)
Y <- as(Y, "CsparseMatrix")

cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)
K_vec <- 1:5

fit_list <- list()

for (cc in cc_vec) {
  
  fit_list[[as.character(cc)]] <- list()
  
  for (K in K_vec) {
    
    set.seed(1)
    fit <- fit_poisson_log1p_nmf(
      Y = Y,
      K = K,
      cc = cc,
      s = rep(1, n),
      loglik = "exact",
      init_method = "random",
      control = list(maxiter = 1000, tol = 1e-8)
    )
    
    fit_list[[as.character(cc)]][[as.character(K)]] <- fit
    
  }
  
}

K_track <- c()
cc_track <- c()
loglik_track <- c()
L_sparsity_track <- c()
F_sparsity_track <- c()


for (cc in cc_vec) {
  
  for (K in K_vec) {
    
    K_track <- c(K_track, K)
    cc_track <- c(cc_track, cc)
    
    fit <- fit_list[[as.character(cc)]][[as.character(K)]]
    
    loglik_track <- c(loglik_track, logLik(fit, Y))
    L_sparsity_track <- c(L_sparsity_track, mean(fit$LL < 1e-10))
    F_sparsity_track <- c(F_sparsity_track, mean(fit$FF < 1e-10))
    
  }
  
}

out_df <- data.frame(
  cc = cc_track,
  K = K_track,
  loglik = loglik_track,
  L_sparsity = L_sparsity_track,
  F_sparsity = F_sparsity_track
)

library(ggplot2)

ggplot(data = out_df, aes(x = cc, y = loglik, group = as.factor(K))) +
  geom_point(aes(color = as.factor(K))) +
  geom_line(aes(color = as.factor(K))) +
  scale_x_log10()

ggplot(data = out_df, aes(x = cc, y = L_sparsity, group = as.factor(K))) +
  geom_point(aes(color = as.factor(K))) +
  geom_line(aes(color = as.factor(K))) +
  scale_x_log10()

ggplot(data = out_df, aes(x = cc, y = F_sparsity, group = as.factor(K))) +
  geom_point(aes(color = as.factor(K))) +
  geom_line(aes(color = as.factor(K))) +
  scale_x_log10()