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

## the fitting/init/trace machinery (ctrl, fit_one, make_rank1_init, init_*,
## chain_traces, run_job) is shared with the MOCA experiment
source("exact_vs_approx_helpers.R")
