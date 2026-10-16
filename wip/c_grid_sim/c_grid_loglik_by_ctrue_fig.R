# Log-likelihood across fitted c, one panel per simulating c.
#
# y is the log-likelihood SHORTFALL from that seed's own best c, i.e.
#   max_over_c_fit(loglik for this seed) - loglik,  so 0 marks the winner.
# Raw log-likelihood cannot be overlaid across seeds: each seed is a different
# dataset, so the curves would be offset by dataset-level constants of a few
# hundred units and the shape -- the thing being compared -- would be
# unreadable. Subtracting each seed's own maximum removes exactly that constant.
#
# Comparing across c_fit within a seed is fair: all eight are Poisson models for
# the same y under the same dominating measure with constants included.
#
# y is log1p because the shortfall is exactly 0 at the winner (which log10
# cannot show) and differences below ~1 log-likelihood unit are meaningless.
#
# Panels share a y scale on purpose: the penalty for being in the wrong REGIME
# (~1000+) dwarfs the penalty within a regime (~1), and a free y would rescale
# every panel and hide that.

source("c_grid_panel_common.R")

b <- c_grid_load()
p <- c_grid_panels(
  b, "subopt",
  ylab     = "log-likelihood shortfall from that seed's best c  (0 = best)",
  title    = "Cost of fitting at the wrong c",
  subtitle = SUB_SEEDS,
  trans = "log1p", breaks = c(0, 1, 10, 100, 1000),
  labels = c("0", "1", "10", "100", "1000"))

ggsave("../images/c_grid_loglik_by_ctrue.png", p, width = 12, height = 6.4,
       dpi = 150, bg = "white")
message("Wrote ../images/c_grid_loglik_by_ctrue.png")
