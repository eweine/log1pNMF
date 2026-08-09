# After the array finishes: which (seed, cc, method) fits are missing, and
# which ran but did not converge? Run from this directory on the cluster:
#
#   Rscript check_approx_sim.R
#
# Prints the missing task ids in sbatch --array form, so a resubmission is
#
#   sbatch --array=<printed list> run_approx_sim.sbatch

source("approx_sim_common.R")

ids     <- 0:(N_TASKS - 1L)
tags    <- vapply(ids, function(i) tag_of(task_spec(i)), character(1))
files   <- file.path(OUT_DIR, paste0(tags, ".rds"))
have    <- file.exists(files)

cat(sum(have), "of", N_TASKS, "outputs present in", OUT_DIR, "\n")

if (any(!have)) {
  miss <- ids[!have]
  cat("\nMISSING task ids (resubmit with sbatch --array=...):\n")
  cat(paste(miss, collapse = ","), "\n")
  for (i in miss) {
    s <- task_spec(i)
    cat(sprintf("  %4d  seed %2d  c_true %-5s  c_fit %-6s  %s\n",
                i, s$seed, fmtc(s$c_true), fmtc(s$cc), s$method))
  }
}

if (any(have)) {
  rows <- do.call(rbind, lapply(files[have], function(f) readRDS(f)$row))
  bad  <- rows[!rows$converged, , drop = FALSE]
  cat("\n", nrow(rows), "fits read;", sum(!rows$converged), "did not converge\n")
  if (nrow(bad) > 0) {
    cat("UNCONVERGED (hit MAXITER):\n")
    print(bad[, c("task", "seed", "c_true", "cc", "method", "n_iter", "seconds")],
          row.names = FALSE)
  }
  cat("\nseconds per fit by method (median [max]):\n")
  for (m in METHODS) {
    sec <- rows$seconds[rows$method == m]
    if (length(sec)) cat(sprintf("  %-10s %8.1f  [%8.1f]\n", m, median(sec), max(sec)))
  }
}
