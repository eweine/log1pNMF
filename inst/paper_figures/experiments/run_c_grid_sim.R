# Array worker for the c-grid simulation: one (seed, c_true, c_fit) cell per
# invocation, single-threaded.
#
#   Rscript run_c_grid_sim.R            # reads SLURM_ARRAY_TASK_ID (0-based)
#   Rscript run_c_grid_sim.R 137        # or take the task id on the command line
#
# Each task regenerates its dataset from (seed, c_true) -- generation is
# deterministic and cheap, so the eight c_fit tasks sharing a dataset build
# byte-identical counts without any shared state or staging step. Each writes
# one small .rds, so a task that dies costs only its own cell.

suppressMessages({
  library(log1pNMF)
  library(fastTopics)
  library(Matrix)
})

source("c_grid_sim_common.R")

## ---- task id ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
id <- if (length(args) >= 1L) as.integer(args[1]) else
        as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(id)) stop("no task id: pass one as an argument or set SLURM_ARRAY_TASK_ID")

spec <- task_spec(id)
tag  <- tag_of(spec)
message(sprintf("=== task %d: %s  (seed %d, c_true %s, c_fit %s, init %s) on %s ===",
                id, tag, spec$seed, format(spec$c_true), format(spec$c_fit),
                spec$init, Sys.info()[["nodename"]]))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT_DIR, paste0(tag, ".rds"))

## ---- data ------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
d  <- sim_dataset(spec$seed, spec$c_true)
message(sprintf("data: %d x %d, %.1f%% zeros, mean %.2f, max %d  (%.1fs)",
                nrow(d$Y), ncol(d$Y), 100 * mean(d$Y == 0), mean(d$Y), max(d$Y),
                proc.time()[["elapsed"]] - t0))

## ---- fit -------------------------------------------------------------------
## c_fit = Inf is the c -> Inf limit of the family, i.e. plain Poisson NMF; there
## is no finite cc that fit_poisson_log1p_nmf can represent it with, so that cell
## is fit with fastTopics. Everything else uses the exact log1p objective.
set.seed(1000L * spec$seed + id)
t1 <- proc.time()[["elapsed"]]

## Both fitters run to convergence on the same criterion -- an absolute
## log-likelihood change below TOL between successive iterations -- with MAXITER
## and MAX_TIME as safety caps only. fit_cell() cold-starts every c_fit,
## including Inf, from a rank-1 warm start; see c_grid_sim_common.R.
f <- fit_cell(d$Y, spec$c_fit)
LL_fit <- f$LL; FF_fit <- f$FF; lam_fit <- f$lam
n_iter <- f$n_iter; converged <- f$converged; fitter <- f$fitter

elapsed <- proc.time()[["elapsed"]] - t1
message(sprintf("fit: %s, %d iterations, converged=%s  (%.1fs)",
                fitter, n_iter, converged, elapsed))

## ---- score and save --------------------------------------------------------
score <- score_fit(d, LL_fit, FF_fit, lam_fit)

row <- cbind(
  data.frame(task = id, seed = spec$seed, c_true = spec$c_true,
             c_fit = spec$c_fit, init = spec$init,
             K_true = K_TRUE, K_fit = K_FIT,
             n = N_ROWS, p = N_COLS,
             alpha = d$alpha, beta = d$beta,
             zero_frac = mean(d$Y == 0), mean_Y = mean(d$Y), max_Y = max(d$Y),
             fitter = fitter, n_iter = n_iter, converged = converged,
             seconds = elapsed, stringsAsFactors = FALSE),
  score)

saveRDS(list(row = row, LL = LL_fit, FF = FF_fit,
             L_true = d$truth$L, F0_true = d$truth$F0, FF_true = d$FF_true),
        out_file)

message("Wrote ", out_file)
print(row[, c("c_true", "c_fit", "init", "L_cor_mean", "F_cor_mean",
              "rate_kl", "loglik", "loglik_oracle", "n_iter", "seconds")])
