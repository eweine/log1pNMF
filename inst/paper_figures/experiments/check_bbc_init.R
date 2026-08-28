# Coverage / convergence check for the BBC initialization-stability
# experiment. Run on the server after (or during) the arrays:
#
#   Rscript check_bbc_init.R
#
# Prints which of the 93 tasks are missing or unconverged, a per-(cc, init)
# summary of the exact log-likelihoods, and ready-to-paste resubmission
# commands for whatever is absent. Reading the outputs is cheap (they are
# small), so there is no caching.

source("bbc_init_common.R")

files <- file.path(OUT_DIR, paste0(sapply(seq_len(N_TASKS),
  function(i) tag_of(SPEC[i, ])), ".rds"))
have  <- file.exists(files)
cat(sprintf("%d / %d outputs present in %s\n", sum(have), N_TASKS, OUT_DIR))

if (any(have)) {
  rows <- do.call(rbind, lapply(files[have], function(f) readRDS(f)$row))
  bad  <- rows[!rows$converged, ]
  if (nrow(bad) > 0) {
    cat("\nUNCONVERGED (hit the iteration cap; consider rerunning with a",
        "higher cap or accepting as-is):\n")
    print(bad[, c("task", "cc", "init", "seed", "n_iter", "elapsed_s")])
  }
  cat("\nExact log-likelihood by (cc, init):\n")
  agg <- do.call(rbind, lapply(split(rows, list(rows$cc, rows$init),
                                     drop = TRUE), function(d)
    data.frame(cc = fmtc(d$cc[1]), init = d$init[1], n = nrow(d),
               best = max(d$loglik_exact), worst = min(d$loglik_exact),
               spread = max(d$loglik_exact) - min(d$loglik_exact),
               median_iter = stats::median(d$n_iter),
               median_hours = stats::median(d$elapsed_s) / 3600)))
  print(agg[order(agg$cc, agg$init), ], row.names = FALSE, digits = 6)
}

miss <- which(!have) - 1L                 # 0-based task ids
if (length(miss) > 0) {
  cat("\nMissing task ids:", paste(miss, collapse = ", "), "\n")
  blk <- function(ids) {                  # contiguous ranges for --array
    if (length(ids) == 0) return(character(0))
    br <- cumsum(c(1L, diff(ids) != 1L))
    sapply(split(ids, br), function(g)
      if (length(g) == 1) as.character(g) else paste0(g[1], "-", g[length(g)]))
  }
  small <- miss[miss <= 30];  big <- miss[miss > 30]
  if (length(small) > 0)
    cat(sprintf("  sbatch --array=%s --cpus-per-task=8 run_bbc_init.sbatch\n",
                paste(blk(small), collapse = ",")))
  if (length(big) > 0)
    cat(sprintf("  sbatch --array=%s --cpus-per-task=1 run_bbc_init.sbatch\n",
                paste(blk(big), collapse = ",")))
} else {
  cat("\nAll tasks present.\n")
}
