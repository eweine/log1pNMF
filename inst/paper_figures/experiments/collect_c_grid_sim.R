# Assemble the per-task .rds files from the c-grid simulation into one tidy
# table, and report which cells are missing.
#
#   Rscript collect_c_grid_sim.R [out.rds]
#
# Writes c_grid_sim_results.rds (a data frame, one row per (seed, c_true,
# c_fit, init)) into OUT_DIR. The per-task files keep the fitted LL / FF, so anything
# needing the factors themselves (structure plots, say) should read those.

source("c_grid_sim_common.R")

out <- commandArgs(trailingOnly = TRUE)
out <- if (length(out)) out[1] else file.path(OUT_DIR, "c_grid_sim_results.rds")

specs <- lapply(seq_len(N_TASKS) - 1L, task_spec)
files <- file.path(OUT_DIR, paste0(vapply(specs, tag_of, character(1)), ".rds"))
ok    <- file.exists(files)

if (any(!ok)) {
  message(sprintf("MISSING %d of %d cells:", sum(!ok), length(ok)))
  miss <- which(!ok) - 1L
  print(utils::head(data.frame(task = miss,
    seed   = vapply(specs[!ok], function(s) s$seed,   numeric(1)),
    c_true = vapply(specs[!ok], function(s) s$c_true, numeric(1)),
    c_fit  = vapply(specs[!ok], function(s) s$c_fit,  numeric(1))), 40))
} else {
  message("all ", length(ok), " cells present")
}

res <- do.call(rbind, lapply(files[ok], function(f) readRDS(f)$rows))
rownames(res) <- NULL
saveRDS(res, out)
message("Wrote ", out, "  (", nrow(res), " rows = ", sum(ok), " cells x ", length(INITS), " inits)")

## quick look. Each cell was fit from both a rank-1 and a random start; the
## headline tables use the BETTER of the two by log-likelihood, which is the
## sensible summary for a nonconvex objective. The per-init tables below show
## whether that choice mattered.
if (nrow(res)) {
  key  <- paste(res$seed, res$c_true, res$c_fit, sep = "|")
  best <- res[order(key, -res$loglik), ]
  best <- best[!duplicated(paste(best$seed, best$c_true, best$c_fit, sep = "|")), ]

  cat("\nBEST-OF-BOTH-INITS, mean over seeds (rows = c_true, cols = c_fit)\n")
  cat("\nL_cor_mean:\n")
  print(round(tapply(best$L_cor_mean, list(best$c_true, best$c_fit), mean), 3))
  cat("\nrate_kl:\n")
  print(signif(tapply(best$rate_kl, list(best$c_true, best$c_fit), mean), 3))

  cat("\nwhich init won, by c_fit:\n")
  print(table(c_fit = best$c_fit, init = best$init))

  cat("\nmean L_cor_mean by init (rows = c_true, cols = c_fit):\n")
  for (i in INITS) {
    cat("  init = ", i, "\n", sep = "")
    sub <- res[res$init == i, ]
    print(round(tapply(sub$L_cor_mean, list(sub$c_true, sub$c_fit), mean), 3))
  }

  cat("\nnot converged, by c_fit and init:\n")
  print(table(c_fit = res$c_fit, init = res$init, converged = res$converged))

  cat("\nmedian seconds per fit, by c_fit and init:\n")
  print(round(tapply(res$seconds, list(res$c_fit, res$init), stats::median), 1))
}
