# Array worker for the exact-vs-approx simulation: ONE (seed, cc, method) fit
# per invocation, single-threaded, rank-1 initialization, run to convergence.
#
#   Rscript run_approx_sim.R          # reads SLURM_ARRAY_TASK_ID (0-based)
#   Rscript run_approx_sim.R 42       # or take the task id on the command line
#
# APPROX_SIM_TASK_OFFSET (default 0) is added to the array task id, so the
# 1980-task grid can be submitted as two <=1001-task Slurm arrays.
#
# Each task regenerates its own data from the seed (generation is
# deterministic and cheap) and writes one small .rds, so tasks share no state
# and a failure costs one fit.
#
# c_fit = Inf is plain Poisson NMF, fit with fastTopics (SCD updates,
# min.delta.loglik = TOL), initialized from fastTopics' exact rank-1
# solution padded to K -- the same convention as the c-grid experiment
# (c_grid_sim_common.R), so "rank-1 init, own objective, run to
# convergence" means the same thing at every point of the grid.
#
# APPROX_SIM_MAXITER overrides MAXITER, for local smoke tests only.

source("approx_sim_common.R")

## ---- task id ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
id <- if (length(args) >= 1L) as.integer(args[1]) else
        as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(id)) stop("no task id: pass one as an argument or set SLURM_ARRAY_TASK_ID")
id <- id + as.integer(Sys.getenv("APPROX_SIM_TASK_OFFSET", "0"))

maxiter <- as.integer(Sys.getenv("APPROX_SIM_MAXITER", MAXITER))

spec <- task_spec(id)
tag  <- tag_of(spec)
message(sprintf("=== task %d: %s  (seed %d, c_true %s, c_fit %s, method %s) on %s ===",
                id, tag, spec$seed, fmtc(spec$c_true), fmtc(spec$cc),
                spec$method, Sys.info()[["nodename"]]))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT_DIR, paste0(tag, ".rds"))

## The grid includes the 810 round-one fits (identical generator, grid and
## tags), which are reused, not refit: skip any task whose output exists.
## This check runs BEFORE the libraries load, so a skipped task costs only
## R startup (~2s), not the 10-15s of package loading.
## APPROX_SIM_OVERWRITE=1 forces a refit (for deliberate reruns only).
if (file.exists(out_file) && Sys.getenv("APPROX_SIM_OVERWRITE", "0") != "1") {
  message("output exists, skipping: ", out_file)
  quit(save = "no")
}

suppressMessages({
  library(log1pNMF)
  library(fastTopics)
  library(Matrix)
})

## ---- data ------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
d  <- sim_dataset(spec$seed, spec$c_true)
Y  <- as(d$Y, "CsparseMatrix")
message(sprintf("data: %d x %d, %.1f%% zeros, mean %.3f, max %d  (%.1fs)",
                nrow(Y), ncol(Y), 100 * mean(d$Y == 0), mean(d$Y), max(d$Y),
                proc.time()[["elapsed"]] - t0))

## ---- fit -------------------------------------------------------------------
## Rank-1 initialization ALWAYS: for finite c the package runs its own K = 1
## warm-up under the fit's own objective (init_method = "rank1"), so each
## method pays its own initialization -- the same convention as the
## exact-vs-approx timing experiments on the real datasets. For c = Inf the
## analogue is fastTopics' exact rank-1 Poisson NMF solution padded to K
## with RANK1_PAD, following c_grid_sim_common.R. Seeded so any stochastic
## element of the warm-up is reproducible per (seed, cc, method).
##
## The convergence criterion is an absolute change in the fit's own objective
## below TOL = 1e-6; MAXITER is a high safety cap and there is no
## optimization-time budget.
RANK1_PAD <- 1e-8

set.seed(1000L * spec$seed + match(spec$cc, CC_GRID) * 10L +
         match(spec$method, METHODS))
t1 <- proc.time()[["elapsed"]]

if (is.infinite(spec$cc)) {
  stopifnot(spec$method == "exact")
  maxiter <- min(maxiter, MAXITER_INF)
  r1 <- fastTopics:::fit_pnmf_rank1(Y)
  L0 <- cbind(r1$L, matrix(RANK1_PAD, nrow(Y), K_FIT - 1L))
  F0 <- cbind(r1$F, matrix(RANK1_PAD, ncol(Y), K_FIT - 1L))
  rownames(L0) <- rownames(Y); rownames(F0) <- colnames(Y)
  colnames(L0) <- colnames(F0) <- paste0("k", seq_len(K_FIT))
  ft <- fastTopics::fit_poisson_nmf(
    Y, fit0 = fastTopics::init_poisson_nmf(Y, F = F0, L = L0),
    numiter = maxiter, method = "scd",
    control = list(nc = THREADS, min.delta.loglik = TOL), verbose = "none")
  fit <- list(LL = ft$L, FF = ft$F,
              objective_trace = as.numeric(ft$progress$loglik))
  n_iter    <- nrow(ft$progress)
  converged <- n_iter < maxiter
  lam_fit   <- tcrossprod(ft$L, ft$F)
} else {
  fit_args <- list(
    Y = Y, K = K_FIT, cc = spec$cc, s = FALSE,
    init_method = "rank1",
    control = list(maxiter = maxiter, tol = TOL, max_time = Inf,
                   verbose = FALSE, threads = THREADS))
  if (spec$method == "exact") {
    fit_args$loglik <- "exact"
  } else {
    fit_args$loglik <- "approx"
    fit_args$approx_technique <- spec$method
  }
  fit <- do.call(fit_poisson_log1p_nmf, fit_args)
  n_iter    <- length(fit$objective_trace) - 1L
  converged <- isTRUE(fit$converged)
  lam_fit   <- stats::fitted(fit)               # s == 1, so this is lambda
}
elapsed <- proc.time()[["elapsed"]] - t1

## The comparison number: the EXACT Poisson log-likelihood of the fitted
## rates, whichever objective produced them.
ll_exact  <- pois_loglik(d$Y, lam_fit)

message(sprintf("fit [%s]: %d iterations, converged=%s, exact ll %.2f  (%.1fs)",
                spec$method, n_iter, converged, ll_exact, elapsed))

row <- data.frame(
  task = id, seed = spec$seed, c_true = spec$c_true, cc = spec$cc,
  method = spec$method,
  K = K_FIT, n = N_ROWS, p = N_COLS,
  zero_frac = mean(d$Y == 0), mean_Y = mean(d$Y), max_Y = max(d$Y),
  n_iter = n_iter, converged = converged, seconds = elapsed,
  loglik_exact = ll_exact,
  loglik_per_entry = ll_exact / (N_ROWS * N_COLS),
  loglik_oracle = pois_loglik(d$Y, d$lambda),
  maxiter_used = maxiter,
  stringsAsFactors = FALSE)

saveRDS(list(row = row,
             LL = as.matrix(fit$LL), FF = as.matrix(fit$FF),
             objective_trace = fit$objective_trace,
             L_true = d$L, F_true = d$FF),
        out_file)

message("Wrote ", out_file)
print(row[, c("seed", "c_true", "cc", "method", "n_iter", "converged",
              "seconds", "loglik_per_entry")], row.names = FALSE)
