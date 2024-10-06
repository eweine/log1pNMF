get_marginal_lik <- function(y, pi0, mu) {

  p <- 1 / (1 + mu)
  lik <- pi0 * ifelse(y == 0, 1, 0) + (1 - pi0) * (1 - p) * (p ^ y)
  total_lik <- sum(log(lik))
  return(total_lik)

}

eb_opt_fn <- function(par, y) {

  pi0 <- boot::inv.logit(par[1])
  mu <- exp(par[2])
  -get_marginal_lik(y, pi0, mu)

}

get_eb_opt <- function(y) {

  opt_out <- optim(
    par = rnorm(2),
    fn = eb_opt_fn,
    y = y
  )

  pi0_out <- boot::inv.logit(opt_out$par[1])
  mu_out <- exp(opt_out$par[2])

  return(
    list(
      pi0 = pi0_out,
      mu = mu_out
    )
  )

}

get_pm <- function(y, pi0, mu) {

  nz_pm <- (y + 1) / (mu + 1)
  p <- 1 / (1 + mu)
  post_pi0 <- ifelse(
    y == 0,
    pi0 / (
      pi0 + (1 - pi0) * (1 - p)
    ),
    0
  )
  pm <- (1 - post_pi0) * nz_pm
  return(pm)

}

solve_pois_mean_id_link <- function(y) {

  eb_par <- get_eb_opt(y)
  pm <- get_pm(y, eb_par$pi0, eb_par$mu)
  return(
    list(
      pi0 = eb_par$pi0,
      mu = eb_par$mu,
      pm = pm
    )
  )

}

n <- 3000
pi0 <- 0.5
mu <- 1
lambda <- rexp(n = n, rate = mu)
idx_0 <- which(rbinom(n = n, size = 1, prob = pi0) == 1)
lambda[idx_0] <- 0
y <- rpois(n = n, lambda = lambda)

pm <- solve_pois_mean_id_link(y)
# Now, I need to be able to get the posterior mean



