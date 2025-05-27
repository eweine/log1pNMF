library(passPCA)
set.seed(1)
n <- 1000
p <- 100

X <- matrix(
  data = rexp(n * p, 1),
  nrow = n,
  ncol = p
)

k <- floor(0.8 * length(X))

# sample k random positions among all entries
idx <- sample(seq_along(X), size = k)

# set those positions to zero
X[idx] <- 0

b <- c(rep(0, 80), runif(20, max = 0.5))

cc <- 1

lambda <- cc * (exp(X %*% b) - 1)
y <- rpois(n, lambda)

X_test <- matrix(
  data = rexp(n * p, 1),
  nrow = n,
  ncol = p
)

k <- floor(0.8 * length(X_test))

# sample k random positions among all entries
idx <- sample(seq_along(X_test), size = k)

# set those positions to zero
X_test[idx] <- 0

lambda_test <- cc * (exp(X_test %*% b) - 1)
y_test <- rpois(n, lambda_test)

cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)

out_list <- list()

for (cc in cc_vec) {
  
  out_list[[as.character(cc)]] <- list()
  
  b_fit <- passPCA:::solve_pois_reg_log1p(
    X = cbind(rep(log(cc), n), X),
    y = y[y != 0],
    y_nz_idx = which(y != 0) - 1,
    s = rep(cc, length(which(y != 0))),
    b = c(1, rep(0.1, p)),
    update_indices = 1:p,
    num_iter = 250,
    alpha = 0.001,
    beta = 0.25
  )[,1][2:(p+1)]
  
  out_list[[as.character(cc)]]$b <- b_fit
  
  out_list[[as.character(cc)]]$train_loglik <- sum(
    dpois(
      y, as.vector(cc * (exp(X %*% b_fit) - 1)), log = TRUE
    )
  )
  
  out_list[[as.character(cc)]]$test_loglik <- sum(
    dpois(
      y_test, as.vector(cc * (exp(X_test %*% b_fit) - 1)), log = TRUE
    )
  )
  
  out_list[[as.character(cc)]]$b_sparsity <- mean(b_fit < 1e-15)
  
}

# should look at the (train) log-likelihood,
# sparsity, and potentially test log-likelihood



res_df <- data.frame(
  cc = cc_vec,
  train_loglik = unlist(purrr::map_depth(out_list, 1, function(x){x$train_loglik})),
  test_loglik = unlist(purrr::map_depth(out_list, 1, function(x){x$test_loglik})),
  b_sparsity = unlist(purrr::map_depth(out_list, 1, function(x){x$b_sparsity}))
)

library(ggplot2)

g1 <- ggplot(data = res_df, aes(x = cc, y = train_loglik)) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  ylab("Train Loglik") +
  xlab("c")

g2 <- ggplot(data = res_df, aes(x = cc, y = test_loglik)) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  ylab("Test Loglik") +
  xlab("c")

g3 <- ggplot(data = res_df, aes(x = cc, y = b_sparsity)) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  ylab("Sparsity") +
  xlab("c")

library(ggpubr)

ggarrange(g1, g2, g3, nrow = 2, ncol = 2)

df2 <- data.frame(
  loglik = c(res_df$test_loglik, res_df$train_loglik),
  cc = rep(res_df$cc, 2),
  type = c(rep("test", length(cc_vec)), rep("train", length(cc_vec)))
)


ggplot(data = df2, aes(x = cc, y = loglik)) +
  geom_point(aes(color = type)) +
  geom_line(aes(color = type)) +
  scale_x_log10() +
  xlab("c") +
  ylab("Loglik")
