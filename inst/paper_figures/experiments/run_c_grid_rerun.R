# Array worker for the c-grid rerun: one entry of the worklist per invocation.
#
#   Rscript run_c_grid_rerun.R          # reads SLURM_ARRAY_TASK_ID (0-based)
#   Rscript run_c_grid_rerun.R 12       # or take the worklist index directly
#
# The index is into c_grid_rerun_todo.rds (built by plan_c_grid_rerun.R), NOT
# into the original 0..639 task space -- the original task id travels along in
# the worklist so output lands on the right filename.
#
# mode == "fresh"  refit from scratch with the rank-1 cold start (this is what
#                  fixes the c_fit = Inf column, and what reruns cells that
#                  never produced a file).
# mode == "warm"   resume the saved LL / FF with a fresh iteration and time
#                  budget. Exact for the log1p path; n_iter accumulates across
#                  the two runs.
#
# The previous .rds is moved to OUT_DIR/superseded/ rather than deleted, so a
# rerun is reversible and the old and new fits can be compared.

suppressMessages({
  library(log1pNMF)
  library(fastTopics)
  library(Matrix)
})

source("c_grid_sim_common.R")

args <- commandArgs(trailingOnly = TRUE)
idx <- if (length(args) >= 1L) as.integer(args[1]) else
         as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(idx))
  stop("no worklist index: pass one as an argument or set SLURM_ARRAY_TASK_ID.\n",
       "  (did you submit without --array=0-N ?)")

if (!file.exists(WORKLIST))
  stop("no worklist at ", WORKLIST, " -- run plan_c_grid_rerun.R first")
work <- readRDS(WORKLIST)

if (idx < 0 || idx >= nrow(work))
  stop("index ", idx, " out of range: worklist has ", nrow(work), " entries (0..",
       nrow(work) - 1L, ")")

w    <- work[idx + 1L, ]
spec <- list(seed = w$seed, c_true = w$c_true, c_fit = w$c_fit)
tag  <- tag_of(spec)
out_file <- file.path(OUT_DIR, paste0(tag, ".rds"))

message(sprintf("=== rerun %d/%d: %s  [%s / %s] on %s ===",
                idx, nrow(work) - 1L, tag, w$reason, w$mode,
                Sys.info()[["nodename"]]))

## ---- previous fit (for warm starts and for provenance) ---------------------
prev <- if (file.exists(out_file))
  tryCatch(readRDS(out_file), error = function(e) NULL) else NULL

if (w$mode == "warm" && is.null(prev))
  stop("mode 'warm' but no readable previous fit at ", out_file)

## ---- data ------------------------------------------------------------------
d <- sim_dataset(spec$seed, spec$c_true)
message(sprintf("data: %d x %d, %.1f%% zeros, mean %.2f, max %d",
                nrow(d$Y), ncol(d$Y), 100 * mean(d$Y == 0), mean(d$Y), max(d$Y)))

## ---- fit -------------------------------------------------------------------
set.seed(1000L * spec$seed + w$task)
t1 <- proc.time()[["elapsed"]]

f <- if (w$mode == "warm")
  fit_cell(d$Y, spec$c_fit, init_LL = prev$LL, init_FF = prev$FF) else
  fit_cell(d$Y, spec$c_fit)

elapsed <- proc.time()[["elapsed"]] - t1

n_iter_prev <- if (w$mode == "warm") prev$row$n_iter else 0L
message(sprintf("fit: %s, %d new iterations (%d total), converged=%s  (%.1fs)",
                f$fitter, f$n_iter, f$n_iter + n_iter_prev, f$converged, elapsed))

## ---- score -----------------------------------------------------------------
score <- score_fit(d, f$LL, f$FF, f$lam)

row <- cbind(
  data.frame(task = w$task, seed = spec$seed, c_true = spec$c_true,
             c_fit = spec$c_fit, K_true = K_TRUE, K_fit = K_FIT,
             n = N_ROWS, p = N_COLS,
             alpha = d$alpha, beta = d$beta,
             zero_frac = mean(d$Y == 0), mean_Y = mean(d$Y), max_Y = max(d$Y),
             fitter = f$fitter, n_iter = f$n_iter + n_iter_prev,
             converged = f$converged, seconds = elapsed,
             stringsAsFactors = FALSE),
  score,
  ## provenance: what this row replaced, so the rerun is auditable
  data.frame(rerun_reason  = w$reason,
             rerun_mode    = w$mode,
             n_iter_prev   = n_iter_prev,
             loglik_prev   = if (is.null(prev)) NA_real_ else prev$row$loglik,
             seconds_prev  = if (is.null(prev)) NA_real_ else prev$row$seconds,
             stringsAsFactors = FALSE))

if (!is.null(prev))
  message(sprintf("loglik: %.3f -> %.3f  (gain %.3f)",
                  prev$row$loglik, row$loglik, row$loglik - prev$row$loglik))

## ---- save, keeping the superseded fit --------------------------------------
if (!is.null(prev) || file.exists(out_file)) {
  dir.create(SUPERSEDED, showWarnings = FALSE, recursive = TRUE)
  invisible(file.rename(out_file, file.path(SUPERSEDED, paste0(tag, ".rds"))))
}

saveRDS(list(row = row, LL = f$LL, FF = f$FF,
             L_true = d$truth$L, F0_true = d$truth$F0, FF_true = d$FF_true),
        out_file)

message("Wrote ", out_file)
print(row[, c("c_true", "c_fit", "rerun_mode", "L_cor_mean", "F_cor_mean",
              "rate_kl", "loglik", "n_iter", "converged", "seconds")])
