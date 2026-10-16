# Figure 3: accuracy of the sparse computational approximation, from the
# Simulation 2 experiment (approx_sim). For each c_true in {1e-3, 1, Inf}
# and each finite fitting c, the per-entry likelihood ratio of each
# approximate fit (Taylor / Chebyshev) relative to the exact fit on the
# SAME dataset and fitting c: exp{(ll_approx - ll_exact) / np}, where both
# log-likelihoods are the EXACT Poisson log-likelihood of the fitted rates.
# Points are means over replicates (30 seeds at c_true = 1 and Inf; 100 at
# c_true = 1e-3), with 95% t-based confidence intervals.
#
# Formatting follows the original single-run version of this figure
# (theme_cowplot, default two-color hue scale, one panel per c_true,
# dashed reference at 1, shared right legend, A/B/C labels).
#
# Reads res_list$approx_sim_metrics from ../data/experiment_results.Rdata
# (produced by analyze_approx_sim.R, collected by utils/collect_data.R).

library(ggplot2)

load("../data/experiment_results.Rdata")
res <- res_list$approx_sim_metrics
NP  <- unique(res$n * res$p)
stopifnot(length(NP) == 1)

## per-entry LR of each approximate fit vs the exact fit at the same
## (seed, c_true, c_fit)
key <- with(res, paste(seed, c_true, cc))
exi <- res$method == "exact"
res$lr_pe <- exp((res$loglik_exact -
                  res$loglik_exact[exi][match(key, key[exi])]) / NP)
ap <- subset(res, method != "exact")
ap$approx_type <- ifelse(ap$method == "taylor", "Taylor", "Chebyshev")

## mean and 95% t-based CI per (c_true, c_fit, method)
agg <- do.call(rbind, lapply(
  split(ap, list(ap$c_true, ap$cc, ap$approx_type), drop = TRUE),
  function(d) data.frame(
    c_true = d$c_true[1], cc = d$cc[1], approx_type = d$approx_type[1],
    n = nrow(d), lr = mean(d$lr_pe),
    ci = stats::qt(0.975, nrow(d) - 1) * stats::sd(d$lr_pe) / sqrt(nrow(d)))))

fmt_c <- function(x) ifelse(is.infinite(x), "∞",
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))

panel_of <- function(ct) {
  d <- subset(agg, c_true == ct)
  ggplot(d, aes(x = cc, y = lr, color = approx_type)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_errorbar(aes(ymin = lr - ci, ymax = lr + ci), width = 0.25) +
    geom_point() +
    geom_line() +
    cowplot::theme_cowplot() +
    xlab("c Used for Fitting") +
    ylab("Per-Entry Likelihood Ratio") +
    ggtitle(paste0("Data Simulated with c = ", fmt_c(ct))) +
    theme(
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 14.5)
    ) +
    ylim(0, 1.05) +
    scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4),
                       transform = "log10") +
    labs(color = "Approximation")
}

g1 <- panel_of(1e-3)
g2 <- panel_of(1)
g3 <- panel_of(Inf)

library(ggpubr)
g <- ggarrange(g1, g2, g3, nrow = 1, common.legend = TRUE,
               legend = "right", labels = "AUTO")

ggsave(
  "../images/approx_quality.png",
  g,
  device = "png",
  width = 12.25,
  height = 3.5,
  bg = "white"
)
message("Wrote ../images/approx_quality.png")
