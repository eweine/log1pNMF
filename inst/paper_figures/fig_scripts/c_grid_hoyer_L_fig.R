# Mean Hoyer sparsity of the fitted loadings, one panel per simulating c.
#
# Hoyer sparsity of a length-n column is (sqrt(n) - L1/L2) / (sqrt(n) - 1): 0 when
# every entry is equal, 1 when a single entry carries everything. Averaged over
# the K = 3 columns of L. The dashed line is the same quantity on the TRUE L,
# identical in every panel because L is drawn once per seed and shared across
# c_true.
#
# Written once per initialization. The two inits agree to within 0.008
# log-likelihood everywhere, but they do NOT produce equivalent fits: 29% of
# cells differ by more than 0.01 in this kind of metric, and at c_fit = 1e-3
# rank-1 gives about twice the factor sparsity of random. Picking by
# log-likelihood would therefore report the random answer at the small-c end
# (it wins all 80 ties there) and a mixture elsewhere, so the init is fixed
# explicitly instead.

source("c_grid_panel_common.R")

for (ini in names(INIT_LAB)) {
  b <- c_grid_load(init = ini)
  p <- c_grid_panels(
    b, "hoyer_L", ref = "hoyer_L_true",
    ylab     = "mean Hoyer sparsity of loadings",
    title    = paste0("Sparsity of the recovered loadings \u2014 ", INIT_LAB[[ini]]),
    subtitle = SUB_REF)
  f <- sprintf("../images/c_grid_hoyer_L_%s.png", ini)
  ggsave(f, p, width = 12, height = 6.4, dpi = 150, bg = "white")
  message("Wrote ", f)
}
