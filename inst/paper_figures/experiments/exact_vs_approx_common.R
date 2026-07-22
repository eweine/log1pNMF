# Shared setup for the exact-vs-approx initialization-timing experiment.
#
# Sourced by the three scheme scripts (exact_vs_approx_scheme{1,2,3}_*.R), each of
# which can be launched independently so all three run in parallel on a server.
# Keeping the config here guarantees the parallel runs stay mutually consistent.
#
# The experiment fits the log1p Poisson NMF model twice from an *identical*
# starting point -- once optimizing the approximate objective and once the exact
# objective -- and records, at every iteration, the cumulative wall-clock time and
# the exact Poisson log-likelihood (via track_time / track_exact_loglik). The
# three schemes differ only in how that shared starting point is built.

library(Matrix)
library(log1pNMF)

## ---- configuration (keep identical across the parallel scheme runs) ---------
## The data file is expected to load a sparse counts matrix (rows = observations,
## cols = features) into the variable named by `data_var`; matches fit_lsa_models.R.
data_path    <- "/rafalab/eweine/log1p_experiments/pancreas_cytokine_lsa.Rdata"  # loads `counts`
data_var     <- "counts"
K            <- 13          # rank of the factorization
cc           <- 1           # link-function tuning parameter c
maxiter      <- 250         # iterations for each full fit
init_maxiter <- 5           # rank-1 warm-up iterations (package default)
seed         <- 1
n_threads    <- 48
out_dir      <- "/rafalab/eweine/log1p_experiments/"    # where the .rds outputs go
out_tag      <- sprintf("exact_vs_approx_init_timing_K%d_c%s", K, format(cc))
tol          <- 1e-8

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
  maxiter = maxiter, init_maxiter = init_maxiter, seed = seed,
  n_threads = n_threads, tol = tol, n = n, p = p
)

## verbose = TRUE so the per-iteration progress bar is visible in the server log
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

## Fit both methods from a shared, explicitly provided initialization.
fit_both <- function(init_LL, init_FF, label) {
  message("  [", label, "] fitting approximate objective ...")
  approx <- fit_poisson_log1p_nmf(
    Y = Y, cc = cc, loglik = "approx",
    init_LL = init_LL, init_FF = init_FF, control = ctrl(maxiter)
  )
  message("  [", label, "] fitting exact objective ...")
  exact <- fit_poisson_log1p_nmf(
    Y = Y, cc = cc, loglik = "exact",
    init_LL = init_LL, init_FF = init_FF, control = ctrl(maxiter)
  )
  list(approx = approx, exact = exact)
}

## Reproduce the package's rank-1 warm-up as a standalone tracked K = 1 exact fit
## (random rank-1 draw + init_maxiter exact updates), padded to rank K. Recording
## it explicitly captures the warm-up's time / log-likelihood trajectory. It is
## seeded, so schemes 2 and 3 build an identical warm-up in their own runs. (For
## cc = 1 this reproduces the internal rank1 exactly; for cc > 1 it differs only
## by the post-fit sqrt(alpha) scaling, which does not affect the comparison since
## both methods share this init.)
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

## Write one scheme's fitted objects and its tidy from-scratch trace table.
## `prefix` is the list of shared warm-up fits prepended to each method's curve.
save_scheme <- function(scheme, res, prefix = list()) {
  fits_path   <- file.path(out_dir, paste0(out_tag, "_", scheme, "_fits.rds"))
  traces_path <- file.path(out_dir, paste0(out_tag, "_", scheme, "_traces.rds"))
  saveRDS(list(config = config, scheme = scheme, results = res, prefix = prefix),
          fits_path)
  traces <- do.call(rbind, lapply(c("approx", "exact"), function(m) {
    curve <- chain_traces(c(prefix, list(res[[m]])))
    curve$scheme <- scheme
    curve$method <- m
    curve$iter   <- seq_len(nrow(curve)) - 1L
    curve
  }))
  rownames(traces) <- NULL
  saveRDS(traces, traces_path)
  message("Wrote ", fits_path)
  message("Wrote ", traces_path)
}
