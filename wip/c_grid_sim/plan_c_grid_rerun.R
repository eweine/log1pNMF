# Work out which c-grid cells need to be redone, and write a worklist for
# run_c_grid_rerun.R.
#
#   Rscript plan_c_grid_rerun.R
#
# Run this on the cluster (it reads OUT_DIR). It classifies all 640 cells and
# writes c_grid_rerun_todo.rds, then prints the exact sbatch command with the
# right array range.
#
# Two reasons a cell is redone:
#
#   missing      no .rds, or one that will not readRDS (a task that never ran,
#                or was killed mid-write -- e.g. a node that accepted the job
#                and then hung). Refit cold. See blame_nodes.sh for finding
#                which nodes those were.
#   unconverged  ran, but EITHER init reported converged == FALSE. Both are
#                resumed from their own saved LL/FF
#                with a fresh iteration budget, which is exact for the log1p
#                path (see fit_cell in c_grid_sim_common.R), so n_iter
#                accumulates rather than restarting.
#
# "fresh" cells discard the old fit; "warm" cells continue it.

source("c_grid_sim_common.R")

specs <- lapply(seq_len(N_TASKS) - 1L, task_spec)
tags  <- vapply(specs, tag_of, character(1))
files <- file.path(OUT_DIR, paste0(tags, ".rds"))

message("scanning ", OUT_DIR, " ...")

status <- vapply(seq_along(files), function(i) {
  if (!file.exists(files[i])) return("missing")
  r <- tryCatch(readRDS(files[i])$rows, error = function(e) NULL)
  if (is.null(r) || !is.data.frame(r)) return("missing")     # truncated / corrupt
  if (nrow(r) != length(INITS))        return("missing")     # one init never finished
  if (!all(r$converged))               return("unconverged") # either init short
  "ok"
}, character(1))

cat("\ncell status over all ", N_TASKS, " cells:\n", sep = "")
print(table(status))

## where the unconverged ones sit, so it is obvious whether they cluster
if (any(status == "unconverged")) {
  cat("\nunconverged cells (rows = c_true, cols = c_fit):\n")
  u  <- status == "unconverged"
  ct <- vapply(specs, function(s) s$c_true, numeric(1))
  cf <- vapply(specs, function(s) s$c_fit,  numeric(1))
  print(table(c_true = ct[u], c_fit = cf[u]))
}

todo <- which(status != "ok")
if (!length(todo)) {
  message("\nnothing to redo.")
  quit(save = "no")
}

work <- data.frame(
  task   = todo - 1L,
  seed   = vapply(specs[todo], function(s) s$seed,   numeric(1)),
  c_true = vapply(specs[todo], function(s) s$c_true, numeric(1)),
  c_fit  = vapply(specs[todo], function(s) s$c_fit,  numeric(1)),
  reason = status[todo],
  mode   = ifelse(status[todo] == "unconverged", "warm", "fresh"),
  stringsAsFactors = FALSE)

saveRDS(work, WORKLIST)
message("\nWrote ", WORKLIST, "  (", nrow(work), " cells)")
print(table(work$reason, work$mode))

cat("\nsubmit with:\n\n  sbatch --array=0-", nrow(work) - 1L,
    " run_c_grid_rerun.sbatch\n\n", sep = "")
