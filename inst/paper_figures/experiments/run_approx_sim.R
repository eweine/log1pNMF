# Array worker for the exact-vs-approx simulation: ONE (seed, cc, method) fit
# per invocation, single-threaded, rank-1 initialization, run to convergence.
#
#   Rscript run_approx_sim.R          # reads SLURM_ARRAY_TASK_ID (0-based)
#   Rscript run_approx_sim.R 42       # or take the task id on the command line
#
# Each task regenerates its own data from the seed (generation is
# deterministic and cheap) and writes one small .rds, so tasks share no state
# and a failure costs one fit.
#
# APPROX_SIM_MAXITER overrides MAXITER, for local smoke tests only.

suppressMessages({
  library(log1pNMF)
  library(Matrix)
})

source("approx_sim_common.R")

## ---- task id ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
id <- if (length(args) >= 1L) as.integer(args[1]) else
        as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(id)) stop("no task id: pass one as an argument or set SLURM_ARRAY_TASK_ID")

maxiter <- as.integer(Sys.getenv("APPROX_SIM_MAXITER", MAXITER))

spec <- task_spec(id)
tag  <- tag_of(spec)
message(sprintf("=== task %d: %s  (seed %d, c_true %s, c_fit %s, method %s) on %s ===",
                id, tag, spec$seed, fmtc(spec$c_true), fmtc(spec$cc),
                spec$method, Sys.info()[["nodename"]]))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT_DIR, paste0(tag, ".rds"))

## ---- data ------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
d  <- sim_dataset(spec$seed, spec$c_true)
Y  <- as(d$Y, "CsparseMatrix")
message(sprintf("data: %d x %d, %.1f%% zeros, mean %.3f, max %d  (%.1fs)",
                nrow(Y), ncol(Y), 100 * mean(d$Y == 0), mean(d$Y), max(d$Y),
                proc.time()[["elapsed"]] - t0))

## ---- fit -------------------------------------------------------------------
## Rank-1 initialization ALWAYS (init_method = "rank1"): the package runs its
## own K = 1 warm-up under the fit's own objective, so each method pays its
## own initialization -- the same convention as the exact-vs-approx timing
## experiments on the real datasets. Seeded so the warm-up draw is
## reproducible per (seed, cc, method).
##
## The convergence criterion is an absolute change in the fit's own objective
## below TOL = 1e-6; MAXITER is a high safety cap and there is no
## optimization-time budget.
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

set.seed(1000L * spec$seed + match(spec$cc, CC_GRID) * 10L +
         match(spec$method, METHODS))
t1  <- proc.time()[["elapsed"]]
fit <- do.call(fit_poisson_log1p_nmf, fit_args)
elapsed <- proc.time()[["elapsed"]] - t1

n_iter    <- length(fit$objective_trace) - 1L
converged <- isTRUE(fit$converged)
lam_fit   <- stats::fitted(fit)                 # s == 1, so this is lambda

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
             LL = fit$LL, FF = fit$FF,
             objective_trace = fit$objective_trace,
             L_true = d$L, F_true = d$FF),
        out_file)

message("Wrote ", out_file)
print(row[, c("seed", "c_true", "cc", "method", "n_iter", "converged",
              "seconds", "loglik_per_entry")], row.names = FALSE)
