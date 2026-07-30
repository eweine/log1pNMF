# Assemble the per-task .rds files from the c-grid simulation into one tidy
# table, and report which cells are missing.
#
#   Rscript collect_c_grid_sim.R [out.rds]
#
# Writes c_grid_sim_results.rds (a data frame, one row per (seed, c_true,
# c_fit)) into OUT_DIR. The per-task files keep the fitted LL / FF, so anything
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

res <- do.call(rbind, lapply(files[ok], function(f) readRDS(f)$row))
rownames(res) <- NULL
saveRDS(res, out)
message("Wrote ", out, "  (", nrow(res), " rows)")

## quick look: mean loading recovery over seeds, c_true down, c_fit across
if (nrow(res)) {
  cat("\nmean L_cor_mean over seeds (rows = c_true, cols = c_fit):\n")
  m <- tapply(res$L_cor_mean, list(res$c_true, res$c_fit), mean)
  print(round(m, 3))
  cat("\nmean rate_kl over seeds (rows = c_true, cols = c_fit):\n")
  print(signif(tapply(res$rate_kl, list(res$c_true, res$c_fit), mean), 3))
  cat("\nmedian seconds per fit, by c_fit:\n")
  print(round(tapply(res$seconds, res$c_fit, stats::median), 1))
}
