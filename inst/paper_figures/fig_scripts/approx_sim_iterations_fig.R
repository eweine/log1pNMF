# Supplementary figure: iterations to convergence in the approximation
# experiment (Simulation 2 datasets), for the exact objective and the two
# approximate objectives, across the finite fitting-c grid. One panel per
# c_true, styled to match the approximation-quality figure (Figure 3).
#
# The c = Inf fits are excluded: they use a different optimizer
# (fastTopics SCD) whose iteration counts are not comparable to the
# log1p fits' update-cycle counts.
#
# Points: geometric mean over replicates (30 per regime, 100 at
# c_true = 1e-3) with 95% t-based confidence intervals computed on the
# log scale (iteration counts are heavy-tailed and the y-axis is
# logarithmic, so the multiplicative summary is the natural one).
# Reads approx_sim_metrics.rds (analyze_approx_sim.R).

library(ggplot2)

res <- readRDS("../experiments/approx_sim_metrics.rds")
res <- subset(res, is.finite(cc))
res$Optimization <- factor(
  c(exact = "Exact", chebyshev = "Chebyshev", taylor = "Taylor")[res$method],
  levels = c("Exact", "Chebyshev", "Taylor"))

## Exact in black; Chebyshev/Taylor keep their Figure 3 colors (the
## ggplot2 default two-color hue pair)
MCOL <- c("Exact" = "black", "Chebyshev" = "#F8766D", "Taylor" = "#00BFC4")

fmt_c <- function(x) ifelse(is.infinite(x), "∞",
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))

panel_of <- function(ct) {
  d <- subset(res, c_true == ct)
  agg <- do.call(rbind, lapply(split(d, list(d$cc, d$Optimization),
                                     drop = TRUE), function(g)
    { l <- log(g$n_iter)
      m <- mean(l)
      h <- stats::qt(0.975, nrow(g) - 1) * stats::sd(l) / sqrt(nrow(g))
      data.frame(cc = g$cc[1], Optimization = g$Optimization[1],
                 mean = exp(m), lo = exp(m - h), hi = exp(m + h)) }))
  ggplot(agg, aes(cc, mean, colour = Optimization)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
    geom_point(size = 1.6) +
    geom_line(linewidth = 0.8) +
    cowplot::theme_cowplot() +
    scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4),
                       transform = "log10") +
    scale_y_log10(labels = scales::label_comma()) +
    scale_colour_manual(values = MCOL, name = "Optimization") +
    xlab("c Used for Fitting") +
    ylab("Iterations to Convergence") +
    ggtitle(paste0("Data Simulated with c = ", fmt_c(ct))) +
    theme(
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 13, hjust = 0.5)
    )
}

g1 <- panel_of(1e-3)
g2 <- panel_of(1)
g3 <- panel_of(Inf)

library(ggpubr)
g <- ggarrange(g1, g2, g3, nrow = 1, common.legend = TRUE,
               legend = "right", labels = "AUTO")

ggsave(
  "../images/approx_sim_iterations.png",
  g,
  device = "png",
  width = 12.25,
  height = 3.5,
  bg = "white"
)
message("Wrote ../images/approx_sim_iterations.png")
