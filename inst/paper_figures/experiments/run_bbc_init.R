# Array worker for the BBC initialization-stability experiment: ONE
# (cc, init, seed) fit per invocation, run to convergence.
#
#   Rscript run_bbc_init.R          # reads SLURM_ARRAY_TASK_ID (0-based)
#   Rscript run_bbc_init.R 42       # or take the task id on the command line
#
# Threads come from SLURM_CPUS_PER_TASK (default 1): submit the c = 0.001
# block (ids 0-30) with --cpus-per-task=8 and the rest (ids 31-92) with
# --cpus-per-task=1; see run_bbc_init.sbatch.
#
# Finite c fits use log1pNMF with the EXACT log-likelihood and the paper's
# size factors; c = Inf is plain Poisson NMF via fastTopics (SCD,
# min.delta.loglik = TOL). "random" inits are seeded per task; "rank1" is
# the deterministic initialization used for the paper's fits.
#
# BBC_INIT_MAXITER overrides MAXITER, for local smoke tests only.

source("bbc_init_common.R")

## ---- task id ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
id <- if (length(args) >= 1L) as.integer(args[1]) else
        as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(id)) stop("no task id: pass one as an argument or set SLURM_ARRAY_TASK_ID")

maxiter <- as.integer(Sys.getenv("BBC_INIT_MAXITER", MAXITER))
THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

spec <- task_spec(id)
tag  <- tag_of(spec)
message(sprintf("=== task %d: %s  (cc %s, init %s, seed %s, %d threads) on %s ===",
                id, tag, fmtc(spec$cc), spec$init,
                ifelse(is.na(spec$seed), "-", spec$seed), THREADS,
                Sys.info()[["nodename"]]))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT_DIR, paste0(tag, ".rds"))

## Skip-if-exists BEFORE the libraries and the DTM build, so a completed
## task costs ~2s of R startup. BBC_INIT_OVERWRITE=1 forces a refit.
if (file.exists(out_file) && Sys.getenv("BBC_INIT_OVERWRITE", "0") != "1") {
  message("output exists, skipping: ", out_file)
  quit(save = "no")
}

suppressMessages({
  library(log1pNMF)
  library(fastTopics)
  library(Matrix)
})

## ---- data ------------------------------------------------------------------
t0  <- proc.time()[["elapsed"]]
dat <- build_bbc()
Y   <- dat$Y
s   <- dat$s
message(sprintf("DTM built in %.1fs", proc.time()[["elapsed"]] - t0))

## ---- fit -------------------------------------------------------------------
## Seed convention: random inits use set.seed(seed) so each task draws its
## own starting point; the rank-1 init is deterministic but gets set.seed(1)
## anyway (matching the paper's fits) in case any step is stochastic.
set.seed(if (spec$init == "random") spec$seed else 1L)
t1 <- proc.time()[["elapsed"]]

if (is.infinite(spec$cc)) {
  maxiter <- min(maxiter, MAXITER_INF)
  if (spec$init == "rank1") {
    r1 <- fastTopics:::fit_pnmf_rank1(Y)
    L0 <- cbind(r1$L, matrix(RANK1_PAD, nrow(Y), K_FIT - 1L))
    F0 <- cbind(r1$F, matrix(RANK1_PAD, ncol(Y), K_FIT - 1L))
    rownames(L0) <- rownames(Y); rownames(F0) <- colnames(Y)
    colnames(L0) <- colnames(F0) <- paste0("k", seq_len(K_FIT))
    fit0 <- fastTopics::init_poisson_nmf(Y, F = F0, L = L0)
  } else {
    fit0 <- fastTopics::init_poisson_nmf(Y, k = K_FIT,
                                         init.method = "random")
  }
  ft <- fastTopics::fit_poisson_nmf(
    Y, fit0 = fit0, numiter = maxiter, method = "scd",
    control = list(nc = THREADS, min.delta.loglik = TOL), verbose = "none")
  fit <- list(LL = ft$L, FF = ft$F,
              objective_trace = as.numeric(ft$progress$loglik))
  n_iter    <- nrow(ft$progress)
  converged <- n_iter < maxiter
  lam_fit   <- tcrossprod(ft$L, ft$F)
} else {
  fit <- fit_poisson_log1p_nmf(
    Y = Y, K = K_FIT, s = s, cc = spec$cc,
    loglik = "exact", init_method = spec$init,
    control = list(maxiter = maxiter, tol = TOL, max_time = Inf,
                   verbose = FALSE, threads = THREADS))
  n_iter    <- length(fit$objective_trace) - 1L
  converged <- isTRUE(fit$converged)
  lam_fit   <- stats::fitted(fit)   # includes the size factors s
}
elapsed <- proc.time()[["elapsed"]] - t1

## The comparison number: the exact Poisson log-likelihood of the fitted
## rates. (At c = Inf fastTopics' rates exclude no scaling; at finite c
## fitted() already includes s.)
ll_exact <- pois_loglik(Y, lam_fit)

message(sprintf("fit: %d iterations, converged=%s, exact ll %.2f  (%.1fs = %.1fh)",
                n_iter, converged, ll_exact, elapsed, elapsed / 3600))

row <- data.frame(
  task = id, cc = spec$cc, init = spec$init, seed = spec$seed,
  K = K_FIT, n = nrow(Y), p = ncol(Y),
  n_iter = n_iter, converged = converged,
  loglik_exact = ll_exact, elapsed_s = elapsed, threads = THREADS,
  node = Sys.info()[["nodename"]])

saveRDS(list(row = row, LL = fit$LL, FF = fit$FF,
             objective_trace = fit$objective_trace),
        out_file)
message("wrote ", out_file)
