# Timing probe: measure per-iteration wall-clock time of the approximate vs exact
# fits on the MOCA data, so the full-run duration can be extrapolated BEFORE
# launching it. Runs a few iterations of each method from a shared (cheap) random
# init and reports per-iteration time.
#
# Run on a high-memory node (same footprint as the real fit -- see
# exact_vs_approx_moca_common.R; ~96-128 GB). A few EXACT iterations on 1.33M
# cells can still take tens of minutes each, so this probe is not instant.
#
#   Rscript probe_moca_timing.R > moca_probe.log 2>&1

source("exact_vs_approx_moca_common.R")   # loads Y (cells x genes), config, helpers

n_probe <- 5                              # iterations to time per method

## a cheap shared random init (per-iteration *time* does not depend on the init)
set.seed(seed)
init_LL <- matrix(stats::runif(n * K, 1e-8, 0.05), n, K)
init_FF <- matrix(stats::runif(p * K, 1e-8, 0.05), p, K)

probe <- function(loglik) {
  message("\n--- probing ", loglik, " (", n_probe, " iterations) ---")
  el <- system.time(
    fit <- fit_poisson_log1p_nmf(
      Y = Y, cc = cc, loglik = loglik, init_LL = init_LL, init_FF = init_FF,
      control = list(maxiter = n_probe, tol = -Inf, verbose = TRUE,
                     threads = n_threads, track_time = TRUE,
                     track_exact_loglik = TRUE)
    )
  )[["elapsed"]]
  opt_per_iter <- mean(diff(fit$time_trace))   # optimization only (excl. setup + diagnostic)
  cat(sprintf("[%s] total elapsed %.1f s for %d iters  (incl. one-time setup + exact-loglik tracking)\n",
              loglik, el, n_probe))
  cat(sprintf("[%s] optimization-only per iteration: %.1f s (mean of time_trace diffs)\n",
              loglik, opt_per_iter))
  c(elapsed = el, opt_per_iter = opt_per_iter)
}

cat("======== MOCA per-iteration timing probe ========\n")
cat(sprintf("Data: %d cells x %d genes, %s nonzeros, K = %d\n\n",
            n, p, format(length(Y@x), big.mark = ","), K))
a <- probe("approx")
e <- probe("exact")

cat("\n======== extrapolated full-fit time (optimization only, excludes init) ========\n")
cat(sprintf("  approx: %d iters x %.1f s = %.1f h\n",
            maxiter_approx, a[["opt_per_iter"]], maxiter_approx * a[["opt_per_iter"]] / 3600))
cat(sprintf("  exact : %d iters x %.1f s = %.1f h\n",
            maxiter_exact,  e[["opt_per_iter"]], maxiter_exact * e[["opt_per_iter"]] / 3600))
cat("  (add the exact-loglik tracking overhead -- roughly total_elapsed/n_probe minus\n")
cat("   opt_per_iter, per iteration -- plus the rank-1 init cost, for the real wall-clock.)\n")
