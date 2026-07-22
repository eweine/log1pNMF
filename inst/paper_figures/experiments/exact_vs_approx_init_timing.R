# Exact log-likelihood vs. wall-clock time: approximate vs. exact optimization.
#
# For a single data set we fit the log1p Poisson NMF model twice -- once
# optimizing the approximate objective (loglik = "approx") and once optimizing
# the exact objective (loglik = "exact") -- and record, at every iteration, the
# cumulative wall-clock time and the *exact* Poisson log-likelihood (via the
# track_time / track_exact_loglik control options). This lets us compare how
# quickly each method drives the true log-likelihood up in wall-clock terms.
#
# Because the comparison is sensitive to where the two methods start, we insist
# that both methods be initialized from the *identical* point. We consider
# three initialization schemes:
#
#   1. "random"       -- a shared random initialization.
#   2. "rank1"        -- a shared rank-1 warm-up (the package's rank1 procedure),
#                        padded to rank K.
#   3. "rank1_plus1"  -- the rank-1 warm-up plus one iteration of exact fitting
#                        (as in fit_lsa_approximate_supp.R).
#
# For schemes 2 and 3 the warm-up is run once as an explicit, tracked fit so its
# own time / log-likelihood trajectory is recorded and can be prepended to each
# method's curve (both methods share the same warm-up, so it is computed once).
#
# All fitted objects and a tidy trace table are written to disk at the end.

library(Matrix)
library(log1pNMF)

## ---- configuration ---------------------------------------------------------
## Edit these for the target data set / machine. The data file is expected to
## load a sparse counts matrix (rows = observations, cols = features) into the
## variable named by `data_var`; this matches fit_lsa_models.R.
data_path    <- "../data/panc_cyto_lsa.Rdata"  # loads `counts`
data_var     <- "counts"
K            <- 13          # rank of the factorization
cc           <- 1           # link-function tuning parameter c
maxiter      <- 250         # iterations for each full fit
init_maxiter <- 5           # rank-1 warm-up iterations (package default)
seed         <- 1
n_threads    <- max(parallel::detectCores() - 1, 1)
out_dir      <- "."         # where the .rds outputs are written
out_tag      <- sprintf("exact_vs_approx_init_timing_K%d_c%s", K, format(cc))

## Force every fit to run the full `maxiter` iterations so the trajectories are
## complete for plotting (set to the package default 1e-8 to allow early stops).
tol          <- -Inf

## ---- load data -------------------------------------------------------------
message("Loading data from ", data_path)
load(data_path)
Y <- get(data_var)
Y <- as(Y, "CsparseMatrix")
n <- nrow(Y)
p <- ncol(Y)
message(sprintf("Data: %d x %d, %d nonzeros", n, p, length(Y@x)))

ctrl <- function(mi) {
  list(
    maxiter            = mi,
    init_maxiter       = init_maxiter,
    tol                = tol,
    verbose            = FALSE,
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

results <- list()
rank1 <- NULL
onestep <- NULL
fits_path <- file.path(out_dir, paste0(out_tag, "_fits.rds"))

## Checkpoint the raw fitted objects after each scheme so a crash during a long
## run does not lose completed work; the tidy trace table is rebuilt at the end.
checkpoint <- function() {
  saveRDS(
    list(config = config, results = results, rank1 = rank1, onestep = onestep),
    fits_path
  )
}

config <- list(
  data_path = data_path, data_var = data_var, K = K, cc = cc,
  maxiter = maxiter, init_maxiter = init_maxiter, seed = seed,
  n_threads = n_threads, tol = tol, n = n, p = p
)

## ---- scheme 1: shared random initialization --------------------------------
## Drawn exactly as the package's init_method = "random" does (runif in
## [1e-8, 0.05]), but generated here so both methods receive the same matrices.
message("Scheme 1: random initialization")
set.seed(seed)
rand_LL <- matrix(stats::runif(n * K, 1e-8, 0.05), n, K)
rand_FF <- matrix(stats::runif(p * K, 1e-8, 0.05), p, K)
results$random <- fit_both(rand_LL, rand_FF, "random")
checkpoint()

## ---- shared rank-1 warm-up (schemes 2 and 3) -------------------------------
## Reproduce the package's rank1 procedure as a standalone, tracked K = 1 exact
## fit: a random rank-1 draw followed by init_maxiter exact updates. Recording
## it explicitly captures the warm-up's time / log-likelihood trajectory, which
## is otherwise internal. (For cc = 1 this reproduces the internal rank1 exactly;
## for cc > 1 it differs only by the post-fit sqrt(alpha) scaling, which does not
## affect the fairness of the comparison since both methods share this init.)
message("Shared rank-1 warm-up (K = 1, ", init_maxiter, " exact iters)")
set.seed(seed)
rank1 <- fit_poisson_log1p_nmf(
  Y = Y, K = 1, cc = cc, loglik = "exact",
  init_method = "random", control = ctrl(init_maxiter)
)
rank1_LL <- cbind(rank1$LL, matrix(1e-8, n, K - 1))
rank1_FF <- cbind(rank1$FF, matrix(1e-8, p, K - 1))

checkpoint()

## ---- scheme 2: shared rank-1 initialization --------------------------------
message("Scheme 2: rank-1 initialization")
results$rank1 <- fit_both(rank1_LL, rank1_FF, "rank1")
checkpoint()

## ---- scheme 3: shared rank-1 + one exact iteration -------------------------
## One additional exact iteration from the padded rank-1 init (tracked), then
## fit both methods from the result.
message("Shared rank-1 + 1 exact iteration")
onestep <- fit_poisson_log1p_nmf(
  Y = Y, cc = cc, loglik = "exact",
  init_LL = rank1_LL, init_FF = rank1_FF, control = ctrl(1)
)
message("Scheme 3: rank-1 + 1 exact iteration initialization")
results$rank1_plus1 <- fit_both(onestep$LL, onestep$FF, "rank1_plus1")
checkpoint()

## ---- assemble tidy "from-scratch" trace table ------------------------------
## Each fit carries aligned time_trace / exact_loglik_trace of length 1 + niter,
## both starting at the initial iterate (time 0). To build a from-scratch curve
## we chain any shared warm-up phase(s) before the method's own trace, offsetting
## times by the cumulative warm-up time and dropping the duplicated join point.
chain_traces <- function(fits) {
  t_off <- 0
  out <- NULL
  for (f in fits) {
    tt <- f$time_trace + t_off
    ll <- f$exact_loglik_trace
    if (!is.null(out)) {          # drop the leading (duplicate) join point
      tt <- tt[-1]
      ll <- ll[-1]
    }
    out <- rbind(out, data.frame(time = tt, exact_loglik = ll))
    t_off <- t_off + tail(f$time_trace, 1)
  }
  out
}

## prefix warm-up phases shared by both methods within each scheme
scheme_prefix <- list(
  random      = list(),
  rank1       = list(rank1),
  rank1_plus1 = list(rank1, onestep)
)

traces <- do.call(rbind, lapply(names(results), function(sc) {
  do.call(rbind, lapply(c("approx", "exact"), function(m) {
    curve <- chain_traces(c(scheme_prefix[[sc]], list(results[[sc]][[m]])))
    curve$scheme <- sc
    curve$method <- m
    curve$iter   <- seq_len(nrow(curve)) - 1L
    curve
  }))
}))
rownames(traces) <- NULL

## ---- write everything to disk ----------------------------------------------
traces_path <- file.path(out_dir, paste0(out_tag, "_traces.rds"))

message("Writing fits to ", fits_path)
checkpoint()
message("Writing traces to ", traces_path)
saveRDS(traces, traces_path)

message("Done.")
