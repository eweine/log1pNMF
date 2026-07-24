# Shared setup for the exact-vs-approx timing experiment on the pancreas data.
#
# Sourced by run_experiment.R (submitted via run_pancreas.sbatch), once per
# (scheme, method) combination. One fit per invocation means each is checkpointed
# independently: a job hitting the wall-clock limit loses only its own fit.
# Keeping the config here guarantees the parallel runs stay mutually consistent.
#
# The experiment fits the log1p Poisson NMF model under the approx vs exact
# objectives, recording per-iteration wall-clock time and the exact Poisson
# log-likelihood (track_time / track_exact_loglik). The `random` scheme starts
# both methods from one shared, seeded random point; the `rank1` scheme gives
# each method its own objective-appropriate warm-up, whose optimization time is
# chained into the total (so the reported time includes the initialization cost).

library(Matrix)
library(log1pNMF)

## ---- configuration (keep identical across the parallel runs) ----------------
## The data file is expected to load a sparse counts matrix (rows = observations,
## cols = features) into the variable named by `data_var`; matches fit_lsa_models.R.
data_path      <- "/rafalab/eweine/log1pNMF/inst/paper_figures/data/pancreas_cytokine_lsa.Rdata"  # loads `counts`
data_var       <- "counts"
K              <- 13        # rank of the factorization
cc             <- 1         # link-function tuning parameter c
maxiter_approx <- 100000L   # high safety cap; max_time / tol are the real stops
maxiter_exact  <- 100000L   # high safety cap; max_time / tol are the real stops
init_maxiter   <- 5         # rank-1 warm-up iterations (package default)
seed           <- 1
n_threads      <- 64        # keep identical across all datasets for consistency
out_dir        <- "/home/ericweine/log1p_experiments/"    # where the .rds outputs go
out_tag        <- sprintf("exact_vs_approx_pancreas_K%d_c%s", K, format(cc))
tol            <- 1e-8      # a fit that converges within the budget stops early
max_time       <- 10 * 3600 # stop each fit after 10 h of optimization time

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
  n_threads = n_threads, tol = tol, max_time = max_time, n = n, p = p
)

## the fitting/init/trace machinery (ctrl, fit_one, make_rank1_init, init_*,
## chain_traces, run_job) is shared with the MOCA experiment
source("exact_vs_approx_helpers.R")
