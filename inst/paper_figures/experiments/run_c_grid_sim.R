# Array worker for the c-grid simulation: one (seed, c_true, c_fit) cell per
# invocation, single-threaded, fit from BOTH initializations.
#
#   Rscript run_c_grid_sim.R            # reads SLURM_ARRAY_TASK_ID (0-based)
#   Rscript run_c_grid_sim.R 137        # or take the task id on the command line
#
# Both inits live in one task rather than one each. Slurm defaults MaxArraySize
# to 1001, so a 1280-task array is rejected with "Invalid job array
# specification"; 640 is known to work here. It also means the dataset is built
# once instead of twice and both inits provably see identical counts.
#
# Each task regenerates its own data from (seed, c_true) -- generation is
# deterministic and cheap -- and writes one small .rds, so tasks share no state
# and a failure costs one cell.

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
message(sprintf("=== task %d: %s  (seed %d, c_true %s, c_fit %s) on %s ===",
                id, tag, spec$seed, format(spec$c_true), format(spec$c_fit),
                Sys.info()[["nodename"]]))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT_DIR, paste0(tag, ".rds"))

## ---- data ------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
d  <- sim_dataset(spec$seed, spec$c_true)
message(sprintf("data: %d x %d, %.1f%% zeros, mean %.2f, max %d  (%.1fs)",
                nrow(d$Y), ncol(d$Y), 100 * mean(d$Y == 0), mean(d$Y), max(d$Y),
                proc.time()[["elapsed"]] - t0))

## ---- fit from each initialization ------------------------------------------
## Both fitters stop on the same criterion: an absolute log-likelihood change
## below TOL, with MAXITER as the only other cap (no optimization-time budget).
rows <- list(); fits <- list()

for (ini in INITS) {
  set.seed(1000L * spec$seed + 10L * id + match(ini, INITS))
  t1 <- proc.time()[["elapsed"]]
  f  <- fit_cell(d$Y, spec$c_fit, init = ini)
  elapsed <- proc.time()[["elapsed"]] - t1

  message(sprintf("fit [%-6s]: %s, %d iterations, converged=%s  (%.1fs)",
                  ini, f$fitter, f$n_iter, f$converged, elapsed))

  rows[[ini]] <- cbind(
    data.frame(task = id, seed = spec$seed, c_true = spec$c_true,
               c_fit = spec$c_fit, init = ini,
               K_true = K_TRUE, K_fit = K_FIT,
               n = N_ROWS, p = N_COLS,
               alpha = d$alpha, beta = d$beta,
               zero_frac = mean(d$Y == 0), mean_Y = mean(d$Y), max_Y = max(d$Y),
               fitter = f$fitter, n_iter = f$n_iter, converged = f$converged,
               seconds = elapsed, stringsAsFactors = FALSE),
    score_fit(d, f$LL, f$FF, f$lam))

  fits[[ini]] <- list(LL = f$LL, FF = f$FF)
}

rows <- do.call(rbind, rows)
rownames(rows) <- NULL

saveRDS(list(rows = rows, fits = fits,
             L_true = d$truth$L, F0_true = d$truth$F0, FF_true = d$FF_true),
        out_file)

message("Wrote ", out_file)
print(rows[, c("c_true", "c_fit", "init", "L_cor_mean", "F_cor_mean",
               "rate_kl", "loglik", "n_iter", "converged", "seconds")],
      row.names = FALSE)
