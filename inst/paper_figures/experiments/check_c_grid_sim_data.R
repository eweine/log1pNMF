# Local sanity check for the c-grid simulation: generate (but do not fit) the
# datasets for one or more seeds across the whole c grid and print the marginal
# summaries, so the calibration can be tuned before launching the array.
#
#   Rscript check_c_grid_sim_data.R [seed ...]      # default: seed 1
#
# What to look for: zero_frac and q999_lam should be near TARGET_ZERO and
# TARGET_Q for every row (that is what the calibration enforces), while max_Y
# and the gene-level spread should stay within an order of magnitude or so
# across c. If max_Y blows up at the small-c end, lower TARGET_Q or QPROB.

source("c_grid_sim_common.R")

seeds <- as.integer(commandArgs(trailingOnly = TRUE))
if (!length(seeds) || any(is.na(seeds))) seeds <- 1L

for (s in seeds) {
  message("\n===== seed ", s, " =====")
  truth <- sim_truth(s)
  message(sprintf("L: rows sum to 1 (max dev %.1e); mean max-weight per row %.3f",
                  max(abs(rowSums(truth$L) - 1)), mean(apply(truth$L, 1, max))))
  message(sprintf("F0: %.1f%% exact zeros; nonzero magnitudes in [%.3f, %.1f]",
                  100 * mean(truth$F0 == 0),
                  min(truth$F0[truth$F0 > 0]), max(truth$F0)))

  tab <- do.call(rbind, lapply(C_GRID, function(cc) {
    d <- sim_dataset(s, cc, truth)
    data_summary(d)
  }))
  print(format(tab, digits = 3, scientific = TRUE), row.names = FALSE)
}
