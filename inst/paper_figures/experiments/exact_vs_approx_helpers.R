# Shared helper functions for the exact-vs-approx initialization-timing
# experiment. Sourced by a dataset-specific "common" script AFTER that script has
# (a) loaded the packages, (b) defined the config globals -- Y, n, p, K, cc,
# maxiter_approx, maxiter_exact, init_maxiter, seed, n_threads, tol, out_dir,
# out_tag, config -- and (c) loaded the data into `Y` (cells x genes). The
# functions below read those globals at call time.
#
# The experiment fits the log1p Poisson NMF model from an *identical* starting
# point under the approx vs exact objectives, recording per-iteration wall-clock
# time and the exact Poisson log-likelihood (track_time / track_exact_loglik).

## verbose = TRUE so the per-iteration progress bar is visible in the job log
## (e.g. `tail -f` the redirected output while the job runs).
ctrl <- function(mi) {
  list(
    maxiter            = mi,
    init_maxiter       = init_maxiter,
    tol                = tol,
    verbose            = TRUE,
    threads            = n_threads,
    track_time         = TRUE,
    track_exact_loglik = TRUE,
    max_time           = max_time   # optimization-time budget (seconds); Inf = none
  )
}

## Fit a single method from a shared, explicitly provided initialization. The
## approximate and exact objectives run for different iteration budgets.
fit_one <- function(loglik, init_LL, init_FF, label) {
  mi <- if (loglik == "approx") maxiter_approx else maxiter_exact
  message("  [", label, " / ", loglik, "] fitting (", mi, " iterations) ...")
  fit_poisson_log1p_nmf(
    Y = Y, cc = cc, loglik = loglik,
    init_LL = init_LL, init_FF = init_FF, control = ctrl(mi)
  )
}

## Reproduce the package's rank-1 warm-up as a standalone tracked K = 1 fit
## (random rank-1 draw + init_maxiter updates), padded to rank K, and recorded so
## its time / log-likelihood trajectory can be chained into the total runtime.
## The warm-up uses the SAME objective as the main fit (approx or exact), matching
## the package's init_method = "rank1": that is how each method's real
## initialization cost -- which differs between approx and exact -- is counted.
## Seeded, so a given (scheme, method) rebuilds an identical warm-up across runs.
make_rank1_init <- function(method) {
  message("Rank-1 warm-up (K = 1, ", init_maxiter, " ", method, " iters)")
  set.seed(seed)
  rank1 <- fit_poisson_log1p_nmf(
    Y = Y, K = 1, cc = cc, loglik = method,
    init_method = "random", control = ctrl(init_maxiter)
  )
  list(rank1    = rank1,
       rank1_LL = cbind(rank1$LL, matrix(1e-8, n, K - 1)),
       rank1_FF = cbind(rank1$FF, matrix(1e-8, p, K - 1)))
}

## The initialization schemes. Each returns the starting point (LL, FF) and the
## list of warm-up `prefix` fits to prepend to the from-scratch trace.
##  - random: a single shared random start, identical for both methods (seeded),
##    with no warm-up, so the two curves begin from the exact same point.
##  - rank1: a method-appropriate warm-up whose optimization time is chained into
##    the total; the approx and exact curves therefore start from their own
##    (different) warm-up, each having paid its own initialization cost.
init_random <- function() {
  set.seed(seed)
  list(LL = matrix(stats::runif(n * K, 1e-8, 0.05), n, K),
       FF = matrix(stats::runif(p * K, 1e-8, 0.05), p, K),
       prefix = list())
}
init_rank1 <- function(method) {
  r1 <- make_rank1_init(method)
  list(LL = r1$rank1_LL, FF = r1$rank1_FF, prefix = list(r1$rank1))
}

## Chain a fit's aligned time_trace / exact_loglik_trace after any shared warm-up
## phase(s), offsetting times by the cumulative warm-up time and dropping the
## duplicated join point.
chain_traces <- function(fits) {
  t_off <- 0
  out <- NULL
  for (f in fits) {
    tt <- f$time_trace + t_off
    ll <- f$exact_loglik_trace
    if (!is.null(out)) { tt <- tt[-1]; ll <- ll[-1] }
    out <- rbind(out, data.frame(time = tt, exact_loglik = ll))
    t_off <- t_off + tail(f$time_trace, 1)
  }
  out
}

## Write one (scheme, method) job's fitted object and its tidy from-scratch trace
## table. `prefix` is the list of shared warm-up fits prepended to the curve.
run_job <- function(scheme, method) {
  ini <- switch(scheme,
                random = init_random(),
                rank1  = init_rank1(method),
                stop("unknown scheme: ", scheme))
  fit <- fit_one(method, ini$LL, ini$FF, scheme)

  tag <- paste0(out_tag, "_", scheme, "_", method)
  saveRDS(list(config = config, scheme = scheme, method = method,
               fit = fit, prefix = ini$prefix),
          file.path(out_dir, paste0(tag, "_fit.rds")))
  curve <- chain_traces(c(ini$prefix, list(fit)))
  curve$scheme <- scheme
  curve$method <- method
  curve$iter   <- seq_len(nrow(curve)) - 1L
  rownames(curve) <- NULL
  saveRDS(curve, file.path(out_dir, paste0(tag, "_traces.rds")))
  message("Wrote ", tag, "_fit.rds and ", tag, "_traces.rds")
  invisible(NULL)
}
