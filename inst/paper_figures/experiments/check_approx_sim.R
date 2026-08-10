# After the array finishes: which (seed, cc, method) fits are missing, and
# which ran but did not converge? Run from this directory on the cluster:
#
#   Rscript check_approx_sim.R
#
# Prints the missing task ids in sbatch --array form, so a resubmission is
#
#   sbatch --array=<printed list> run_approx_sim.sbatch
#
# Safe to run while the array is still going; "missing" then includes tasks
# that are simply still running. Each output .rds carries the fit's whole
# objective trace, so reading one is slow; the one-row summaries are cached
# in OUT_DIR and repeated checks only read outputs that appeared (or were
# rewritten) since the last run. Set CHECK_CORES > 1 to read new outputs in
# parallel, e.g.
#
#   CHECK_CORES=8 Rscript check_approx_sim.R

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
  have_files <- files[have]
  have_mtime <- as.numeric(file.mtime(have_files))

  ## cache of already-read summary rows, keyed by (file, mtime)
  cache_file <- file.path(OUT_DIR, "approx_sim_rows_cache.rds")
  cache <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
  if (is.null(cache)) {
    ci  <- rep(NA_integer_, length(have_files))
    hit <- rep(FALSE, length(have_files))
  } else {
    ci  <- match(have_files, cache$file)
    hit <- !is.na(ci) & cache$mtime[ci] == have_mtime
  }

  new_files <- have_files[!hit]
  new_rows  <- NULL
  if (length(new_files) > 0) {
    n_cores <- max(1L, as.integer(Sys.getenv("CHECK_CORES", "1")))
    cat(sprintf("\nreading %d new outputs (%d cached, %d cores)...\n",
                length(new_files), sum(hit), n_cores))
    ## a file being saved by a worker right now can fail to read; skip it
    ## this run and it will be picked up (with a fresh mtime) next time
    read_row <- function(f) tryCatch(cbind(readRDS(f)$row, file = f),
                                     error = function(e) NULL)
    t0 <- proc.time()[["elapsed"]]
    chunks <- split(new_files, ceiling(seq_along(new_files) / 25))
    got <- list(); done <- 0L
    for (ch in chunks) {
      rs <- if (n_cores > 1L) parallel::mclapply(ch, read_row, mc.cores = n_cores)
            else lapply(ch, read_row)
      got <- c(got, rs); done <- done + length(ch)
      el <- proc.time()[["elapsed"]] - t0
      cat(sprintf("  %3d / %d  (%.0fs elapsed, ~%.0fs left)\n", done,
                  length(new_files), el,
                  el / done * (length(new_files) - done)))
    }
    unreadable <- sum(vapply(got, is.null, logical(1)))
    if (unreadable > 0)
      cat(unreadable, "output(s) unreadable (probably mid-write); skipped\n")
    new_rows <- do.call(rbind, got)
    if (!is.null(new_rows))
      new_rows$mtime <- have_mtime[!hit][match(new_rows$file, new_files)]
  }

  cache <- rbind(if (any(hit)) cache[ci[hit], , drop = FALSE], new_rows)
  saveRDS(cache, cache_file)
  rows <- cache
  if (is.null(rows)) { cat("no outputs readable yet\n"); quit(save = "no") }

  bad <- rows[!rows$converged, , drop = FALSE]
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
