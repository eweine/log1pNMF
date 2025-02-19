test_that("a basic fit works", {

  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4)
  fit <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    init_method = "random",
    control = list(verbose = FALSE)
  )

  expect_true(all(fit$LL>=0))
  expect_true(all(fit$FF>=0))
  expect_nondecreasing(fit$objective_trace)

})

test_that("big c works", {
  
  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4)
  fit <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    cc = 1000,
    init_method = "random",
    control = list(verbose = FALSE)
  )
  
  expect_true(all(fit$LL>=0))
  expect_true(all(fit$FF>=0))
  expect_nondecreasing(fit$objective_trace)
  
})

test_that("small c works", {
  
  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4)
  fit <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    cc = 1e-4,
    init_method = "random",
    control = list(verbose = FALSE)
  )
  
  expect_true(all(fit$LL>=0))
  expect_true(all(fit$FF>=0))
  expect_nondecreasing(fit$objective_trace)
  
})

test_that("rank1 initialization works", {

  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4)
  fit <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    init_method = "rank1",
    control = list(verbose = FALSE)
  )

  expect_true(all(fit$LL>=0))
  expect_true(all(fit$FF>=0))
  expect_nondecreasing(fit$objective_trace)

})

test_that("exact loglikelihood works", {

  cc <- 0.75
  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4, cc = cc)
  fit <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    cc = cc,
    loglik = "exact",
    control = list(verbose = FALSE)
  )

  expect_true(all(fit$LL>=0))
  expect_true(all(fit$FF>=0))
  expect_nondecreasing(fit$objective_trace)

})

test_that("exact produces better fits than approximate", {

  cc <- 1
  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4, cc = cc)

  fit_exact <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    cc = cc,
    loglik = "exact",
    control = list(verbose = FALSE)
  )
  Lambda_hat <- fitted(fit_exact)

  fit_approx <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    cc = cc,
    loglik = "approx",
    control = list(verbose = FALSE)
  )
  Lambda_hat_approx <- fitted(fit_approx)
  expect_true(all(fit_approx$LL>=0))
  expect_true(all(fit_approx$FF>=0))
  expect_true(all(fit_exact$LL>=0))
  expect_true(all(fit_exact$FF>=0))

  ll_exact <- sum(
    dpois(
      x = as.vector(as.matrix(dat$Y)),
      lambda = as.vector(Lambda_hat),
      log = TRUE
    )
  )

  ll_approx <- sum(
    dpois(
      x = as.vector(as.matrix(dat$Y)),
      lambda = as.vector(Lambda_hat_approx),
      log = TRUE
    )
  )

  testthat::expect_gt(ll_exact, ll_approx)

})

test_that("Chebyshev approximation works", {

  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4)
  fit <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    init_method = "rank1",
    approx_technique = "chebyshev",
    control = list(verbose = FALSE)
  )

  expect_true(all(fit$LL>=0))
  expect_true(all(fit$FF>=0))
  expect_nondecreasing(fit$objective_trace)

})

test_that("true value of c improves fit", {

  cc <- 1
  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4, cc = cc)

  fit_exact <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    cc = cc,
    loglik = "exact",
    control = list(verbose = FALSE)
  )
  Lambda_hat <- fitted(fit_exact)

  fit_misspec <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    cc = 10,
    loglik = "approx",
    control = list(verbose = FALSE)
  )
  Lambda_hat_misspec <- fitted(fit_misspec)

  ll_exact <- sum(
    dpois(
      x = as.vector(as.matrix(dat$Y)),
      lambda = as.vector(Lambda_hat),
      log = TRUE
    )
  )

  ll_misspec <- sum(
    dpois(
      x = as.vector(as.matrix(dat$Y)),
      lambda = as.vector(Lambda_hat_misspec),
      log = TRUE
    )
  )

  testthat::expect_gt(ll_exact, ll_misspec)
  expect_true(all(fit_approx$LL>=0))
  expect_true(all(fit_approx$FF>=0))
  expect_true(all(fit_exact$LL>=0))
  expect_true(all(fit_exact$FF>=0))

})

test_that("providing custom initializations works", {

  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4)
  fit <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    init_method = "random",
    loglik = "exact",
    control = list(verbose = FALSE, maxiter = 5)
  )

  Lambda_hat <- fitted(fit)

  ll_init <- sum(
    dpois(
      x = as.vector(as.matrix(dat$Y)),
      lambda = as.vector(Lambda_hat),
      log = TRUE
    )
  )

  fit2 <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    init_LL = fit$LL,
    init_FF = fit$FF,
    loglik = "exact",
    control = list(verbose = FALSE, maxiter = 5)
  )

  Lambda_hat2 <- fitted(fit2)

  ll_final <- sum(
    dpois(
      x = as.vector(as.matrix(dat$Y)),
      lambda = as.vector(Lambda_hat2),
      log = TRUE
    )
  )

  expect_true(all(fit$LL>=0))
  expect_true(all(fit$FF>=0))
  expect_true(all(fit2$LL>=0))
  expect_true(all(fit2$FF>=0))
  testthat::expect_gt(ll_final, ll_init)

})

test_that("providing custom size factors works", {

  set.seed(1)
  dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4)
  s <- runif(500, 0.9, 1.1)
  fit <- fit_poisson_log1p_nmf(
    Y = dat$Y,
    K = 4,
    s = s,
    init_method = "random",
    control = list(verbose = FALSE)
  )

  expect_nondecreasing(fit$objective_trace)
  expect_equal(fit$s, s)
  expect_true(all(fit$LL>=0))
  expect_true(all(fit$FF>=0))

})
