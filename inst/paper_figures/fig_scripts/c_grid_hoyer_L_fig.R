# Mean Hoyer sparsity of the fitted loadings, one panel per simulating c.
#
# Hoyer sparsity of a length-n column is (sqrt(n) - L1/L2) / (sqrt(n) - 1):
# 0 when every entry is equal, 1 when a single entry carries everything.
# Averaged over the K = 3 columns of L.
#
# The dashed line is the same quantity on the TRUE L. Because L's rows are
# Dirichlet(0.3) draws the truth is fairly peaked, so a curve well below the
# dashed line means the fit spread each cell's weight across factors that the
# truth kept separate.

source("c_grid_panel_common.R")

b <- c_grid_load()
p <- c_grid_panels(
  b, "hoyer_L", ref = "hoyer_L_true",
  ylab     = "mean Hoyer sparsity of loadings",
  title    = "Sparsity of the recovered loadings",
  subtitle = SUB_REF)

ggsave("../images/c_grid_hoyer_L.png", p, width = 12, height = 6.4,
       dpi = 150, bg = "white")
message("Wrote ../images/c_grid_hoyer_L.png")
