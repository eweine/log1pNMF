# Shared setup for the exact-vs-approx initialization-timing experiment.
#
# Sourced by the six job scripts (exact_vs_approx_<scheme>_<method>.R), each of
# which fits ONE method (approx or exact) from ONE initialization scheme and
# writes its own output. Splitting to one fit per job means each fit is
# checkpointed independently: a job hitting the wall-clock limit only loses its
# own fit, and the other five are unaffected. Keeping the config here guarantees
# the parallel runs stay mutually consistent.
#
# The experiment fits the log1p Poisson NMF model from an *identical* starting
# point under approx vs exact objectives, recording per-iteration wall-clock time
# and the exact Poisson log-likelihood (track_time / track_exact_loglik). The
# three schemes differ only in how that shared starting point is built. Because
# the init is seeded, the approx and exact jobs of a scheme build the identical
# starting point in their own runs.

library(Matrix)
library(log1pNMF)

## ---- configuration (keep identical across the parallel runs) ----------------
## The data file is expected to load a sparse counts matrix (rows = observations,
## cols = features) into the variable named by `data_var`; matches fit_lsa_models.R.
data_path      <- "/home/ericweine/log1p_experiments/pancreas_cytokine_lsa.Rdata"  # loads `counts`
data_var       <- "counts"
K              <- 13        # rank of the factorization
cc             <- 1         # link-function tuning parameter c
maxiter_approx <- 300       # iterations for the approximate fit
maxiter_exact  <- 100       # iterations for the exact fit
init_maxiter   <- 5         # rank-1 warm-up iterations (package default)
seed           <- 1
n_threads      <- 48
out_dir        <- "/home/ericweine/log1p_experiments/"    # where the .rds outputs go
out_tag        <- sprintf("exact_vs_approx_init_timing_K%d_c%s", K, format(cc))
tol            <- 1e-8

## ---- load data --------------------------------------------------------------
message("Loading data from ", data_path)
load(data_path)
Y <- get(data_var)
Y <- as(Y, "CsparseMatrix")
n <- nrow(Y)
p <- ncol(Y)
message(sprintf("Data: %d x %d, %d nonzeros", n, p, length(Y@x)))

config <- list(
  data_path = data_path, data_var = data_var, K = K, cc = cc,
  maxiter_approx = maxiter_approx, maxiter_exact = maxiter_exact,
  init_maxiter = init_maxiter, seed = seed,
  n_threads = n_threads, tol = tol, n = n, p = p
)

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
    track_exact_loglik = TRUE
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

## Reproduce the package's rank-1 warm-up as a standalone tracked K = 1 exact fit
## (random rank-1 draw + init_maxiter exact updates), padded to rank K. Recording
## it explicitly captures the warm-up's time / log-likelihood trajectory. It is
## seeded, so the approx and exact jobs of a scheme build an identical warm-up.
make_rank1_init <- function() {
  message("Shared rank-1 warm-up (K = 1, ", init_maxiter, " exact iters)")
  set.seed(seed)
  rank1 <- fit_poisson_log1p_nmf(
    Y = Y, K = 1, cc = cc, loglik = "exact",
    init_method = "random", control = ctrl(init_maxiter)
  )
  list(rank1    = rank1,
       rank1_LL = cbind(rank1$LL, matrix(1e-8, n, K - 1)),
       rank1_FF = cbind(rank1$FF, matrix(1e-8, p, K - 1)))
}

## The three initialization schemes. Each returns the shared starting point
## (LL, FF) and the list of warm-up `prefix` fits to prepend to the from-scratch
## trace. Seeded, so both methods of a scheme get the identical start.
init_random <- function() {
  set.seed(seed)
  list(LL = matrix(stats::runif(n * K, 1e-8, 0.05), n, K),
       FF = matrix(stats::runif(p * K, 1e-8, 0.05), p, K),
       prefix = list())
}
init_rank1 <- function() {
  r1 <- make_rank1_init()
  list(LL = r1$rank1_LL, FF = r1$rank1_FF, prefix = list(r1$rank1))
}
init_rank1_plus1 <- function() {
  r1 <- make_rank1_init()
  message("Shared rank-1 + 1 exact iteration")
  onestep <- fit_poisson_log1p_nmf(
    Y = Y, cc = cc, loglik = "exact",
    init_LL = r1$rank1_LL, init_FF = r1$rank1_FF, control = ctrl(1)
  )
  list(LL = onestep$LL, FF = onestep$FF, prefix = list(r1$rank1, onestep))
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
                random      = init_random(),
                rank1       = init_rank1(),
                rank1_plus1 = init_rank1_plus1(),
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
