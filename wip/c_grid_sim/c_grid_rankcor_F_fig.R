# Mean rank correlation between the fitted factors, one panel per simulating c.
#
# For each fit: Spearman correlation between every pair of the K = 3 factor
# columns, averaged over the 3 distinct pairs. High values mean the recovered
# factors are redundant rather than distinct parts. The dashed line is the same
# quantity on the TRUE F_c; it is identical in every panel because Spearman
# correlation is invariant to the positive affine map F_c = alpha*F0 + beta.
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
    b, "rankcor_F", ref = "rankcor_F_true",
    ylab     = "mean rank correlation between factors",
    title    = paste0("Are the recovered factors distinct? \u2014 ", INIT_LAB[[ini]]),
    subtitle = SUB_REF)
  f <- sprintf("../images/c_grid_rankcor_F_%s.png", ini)
  ggsave(f, p, width = 12, height = 6.4, dpi = 150, bg = "white")
  message("Wrote ", f)
}
