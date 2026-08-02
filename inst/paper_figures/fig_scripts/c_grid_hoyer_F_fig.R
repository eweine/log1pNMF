# Mean Hoyer sparsity of the fitted factors, one panel per simulating c.
#
# Hoyer sparsity of a length-n column is (sqrt(n) - L1/L2) / (sqrt(n) - 1):
# 0 when every entry is equal, 1 when a single entry carries everything.
# Averaged over the K = 3 columns of F.
#
# The dashed reference is the sparsity of the REALIZED truth F_c = alpha*F0 +
# beta, which is what the fit is actually estimating. Note that F_c is barely
# sparse at all -- the beta offset puts a floor under every entry -- so the
# reference sits near 0 even though the underlying F0 is 50% exact zeros. If the
# question is "did the fit recover the sparsity of the generative structure",
# compare against hoyer_F0_true in the metrics table instead; that column holds
# the sparsity of F0 itself.

source("c_grid_panel_common.R")

b <- c_grid_load()
p <- c_grid_panels(
  b, "hoyer_F", ref = "hoyer_F_true",
  ylab     = "mean Hoyer sparsity of factors",
  title    = "Sparsity of the recovered factors",
  subtitle = SUB_REF)

ggsave("../images/c_grid_hoyer_F.png", p, width = 12, height = 6.4,
       dpi = 150, bg = "white")
message("Wrote ../images/c_grid_hoyer_F.png")
