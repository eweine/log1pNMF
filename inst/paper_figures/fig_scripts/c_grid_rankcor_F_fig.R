# Mean rank correlation between the fitted factors, one panel per simulating c.
#
# For each fit: Spearman correlation between every pair of the K = 3 factor
# columns, averaged over the 3 distinct pairs. High values mean the recovered
# factors are redundant rather than distinct parts.
#
# The dashed line is the same quantity computed on the TRUE F_c for that
# simulation, so a curve sitting above it means the fit recovered factors more
# correlated -- less separated -- than the truth actually is.

source("c_grid_panel_common.R")

b <- c_grid_load()
p <- c_grid_panels(
  b, "rankcor_F", ref = "rankcor_F_true",
  ylab     = "mean rank correlation between factors",
  title    = "Are the recovered factors distinct?",
  subtitle = SUB_REF)

ggsave("../images/c_grid_rankcor_F.png", p, width = 12, height = 6.4,
       dpi = 150, bg = "white")
message("Wrote ../images/c_grid_rankcor_F.png")
