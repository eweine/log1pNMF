# I think it could be interesting to compare my ability to
# estimate log1p b here...

#set.seed(1000)

library(flashier)
library(distr)

n <- 1500
p <- 750
K <- 5

l_dist <- UnivarMixingDistribution(
  Dirac(0),
  Exp(3),
  mixCoeff = c(2/3,1/3)
)

f_dist <- UnivarMixingDistribution(
  Dirac(0),
  Exp(3),
  mixCoeff = c(2/3,1/3)
)

l_sampler <- distr::r(l_dist)
f_sampler <- distr::r(f_dist)

get_n_factors_flash <- function(f) {

  n_factors <- f$n_factors

  for (j in 1:n_factors) {

    if (all(f$F_pm[, j] < 1e-12)) {

      n_factors <- n_factors - 1

    }

  }

  return(n_factors)

}

n_sims <- 10

mse_mle_vec <- c()
mse_log1p_old_vec <- c()
mse_log1p_new_mle_vec <- c()
mse_log1p_new_greedy_vec <- c()

factors_log1p_old_vec <- c()
factors_log1p_new_mle_vec <- c()
factors_log1p_new_greedy_vec <- c()

for (i in 1:n_sims) {

  set.seed(i)

  LL <- matrix(
    data = l_sampler(n * K),
    nrow = n,
    ncol = K
  )

  FF <- matrix(
    data = f_sampler(p * K),
    nrow = p,
    ncol = K
  )

  Lambda <- exp(tcrossprod(LL, FF)) - 1
  y <- rpois(n = length(as.vector(Lambda)), lambda = as.vector(Lambda))

  Y <- matrix(
    data = y,
    nrow = n,
    ncol = p
  )

  rownames(Y) <- paste0("cell", 1:n)
  colnames(Y) <- paste0("gene", 1:p)

  B <- tcrossprod(LL, FF)

  rownames(B) <- paste0("cell", 1:n)
  colnames(B) <- paste0("gene", 1:p)

  Y <- Y[rowSums(Y) > 0, ]
  Y <- Y[, colSums(Y) > 0]

  B <- B[rownames(B) %in% rownames(Y), ]
  B <- B[, colnames(B) %in% colnames(Y)]

  Y_trans <- log1p(Y)

  x  <- rpois(1e7, 1/n)
  s1 <- sd(log(x + 1))

  fit_flash_old <- flash(
    data = as(Y_trans, "sparseMatrix"),
    S = s1,
    ebnm_fn = ebnm::ebnm_point_exponential,
    backfit = TRUE
  )

  library(passPCA)

  log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
    Y = as(Y, "sparseMatrix"),
    K = 10,
    approx_range = c(0, 1.25),
    maxiter = 500
  )

  S <- as.matrix(Y / ((1 + Y) ^ 2))
  S <- apply(S, c(1, 2), function(x){if(x == 0){1}else{x}})

  Y_trans2 <- apply(Y_trans, c(1, 2), function(x){if(x == 0){-1}else{x}})

  fit_flash_new_mle <- flash_init(data = Y_trans2, S = sqrt(S))
  fit_flash_new_mle <- flash_factors_init(
    flash = fit_flash_new_mle,
    init = list(
      log1p_fit$U, log1p_fit$V
    ),
    ebnm_fn = ebnm::ebnm_point_exponential
  )

  fit_flash_new_mle <- flash_backfit(fit_flash_new_mle)
  fit_flash_new_greedy <- flash(data = Y_trans2, S = sqrt(S), backfit = TRUE)


  B_mle <- tcrossprod(log1p_fit$U, log1p_fit$V)
  B_log1p_old <- tcrossprod(fit_flash_old$L_pm, fit_flash_old$F_pm)
  B_log1p_new_mle <- tcrossprod(fit_flash_new_mle$L_pm, fit_flash_new_mle$F_pm)
  B_log1p_new_greedy <- tcrossprod(
    fit_flash_new_greedy$L_pm,
    fit_flash_new_greedy$F_pm
  )

  n_factors_log1p_old <- get_n_factors_flash(fit_flash_old)
  n_factors_log1p_new <- get_n_factors_flash(fit_flash_new_mle)
  n_factors_log1p_new_greedy <- get_n_factors_flash(fit_flash_new_greedy)

  mse_mle <- mean((B_mle - B) ^ 2)
  mse_log1p_old <- mean((B_log1p_old - B) ^ 2)
  mse_log1p_new_mle <- mean((B_log1p_new_mle - B) ^ 2)
  mse_log1p_new_greedy <- mean((B_log1p_new_greedy - B) ^ 2)

  mse_mle_vec <- c(mse_mle_vec, mse_mle)
  mse_log1p_old_vec <- c(mse_log1p_old_vec, mse_log1p_old)
  mse_log1p_new_mle_vec <- c(mse_log1p_new_mle_vec, mse_log1p_new_mle)
  mse_log1p_new_greedy_vec <- c(mse_log1p_new_greedy_vec, mse_log1p_new_greedy)

  factors_log1p_old_vec <- c(factors_log1p_old_vec, n_factors_log1p_old)
  factors_log1p_new_mle_vec <- c(factors_log1p_new_mle_vec, n_factors_log1p_new)
  factors_log1p_new_greedy_vec <- c(factors_log1p_new_greedy_vec, n_factors_log1p_new_greedy)

}

# for now let's suppose that this is correct (though I'm not sure it is)
# then I would like to at least objectively evaluate what I have here

factors_res_df <- data.frame(
  factors = c(factors_log1p_old_vec, factors_log1p_new_mle_vec, factors_log1p_new_greedy_vec),
  method = c(rep("flash_old", 10), rep("flash_new_mle_init", 10), rep("flash_new_greedy_init", 10))
)

library(dplyr)
df_summary <- factors_res_df %>%
  group_by(method) %>%
  summarise(
    mean_factors = mean(factors, na.rm = TRUE),
    sd_factors = sd(factors, na.rm = TRUE)
  )

mse_res_df <- data.frame(
  mse = c(mse_mle_vec, mse_log1p_new_greedy_vec, mse_log1p_new_mle_vec, mse_log1p_old_vec),
  method = c(rep("MLE", 10), rep("flash_new_greedy_init", 10), rep("flash_new_mle_init", 10), rep("flash_old", 10))
)

df_summary2 <- mse_res_df %>%
  group_by(method) %>%
  summarise(
    mean_mse = mean(mse, na.rm = TRUE),
    sd_mse = sd(mse, na.rm = TRUE)
  )


