# Appendix figure for Simulation 2, exact fits on seeds 1-30:
#
#   A  RMSE of the fitted rates against the true rates, relative to the
#      same replicate's fit at c_fit = c_true (log scale; dashed: the
#      reference ratio 1)
#   B  mean ABSOLUTE pairwise Spearman correlation between the columns of
#      the fitted F -- the factor-redundancy measure of the real-data
#      analyses (fig:sparsity, same abs convention) -- with the truth's
#      own value dashed
#
# Thin lines: individual seeds; thick line: mean; open circle:
# c_fit = c_true. Reads res_list$approx_sim_metrics (experiment_results.Rdata) (analyze_approx_sim.R).

library(ggplot2)
library(cowplot)

N_SEEDS_SHOW <- 30

load("../data/experiment_results.Rdata")
res <- res_list$approx_sim_metrics
ex  <- subset(res, method == "exact" & seed <= N_SEEDS_SHOW)

fmt_c <- function(x) ifelse(is.infinite(x), "∞",
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
CC  <- sort(unique(ex$cc))
ex$cf <- factor(fmt_c(ex$cc), levels = fmt_c(CC))
CTL <- paste0("Data Simulated with c = ", fmt_c(sort(unique(ex$c_true))))
ex$ct <- factor(paste0("Data Simulated with c = ", fmt_c(ex$c_true)),
                levels = CTL)
ex$is_ref <- (is.finite(ex$c_true) & ex$cc == ex$c_true) |
             (is.infinite(ex$c_true) & is.infinite(ex$cc))

ref <- ex[ex$is_ref, c("seed", "ct", "rmse")]
names(ref)[3] <- "rmse_ref"
ex <- merge(ex, ref, by = c("seed", "ct"))
ex$rel_rmse <- ex$rmse / ex$rmse_ref

SEED_COL <- "#9BB7BD"; MEAN_COL <- "#0E5C6B"; REF_COL <- "#C2410C"

row_panel <- function(v, ylab, ref1 = FALSE, truth_df = NULL, logy = FALSE,
                      strips = FALSE, xaxis = FALSE) {
  d <- ex; d$value <- ex[[v]]
  avg <- aggregate(value ~ cf + ct, d, mean)
  dg  <- merge(unique(d[d$is_ref, c("cf", "ct")]), avg)
  g <- ggplot(d, aes(cf, value))
  if (ref1)
    g <- g + geom_hline(yintercept = 1, linetype = "dashed",
                        colour = REF_COL, linewidth = 0.4)
  if (!is.null(truth_df))
    g <- g + geom_hline(data = truth_df, aes(yintercept = value),
                        linetype = "dashed", colour = REF_COL,
                        linewidth = 0.4)
  g <- g +
    geom_line(aes(group = seed), colour = SEED_COL,
              linewidth = 0.3, alpha = 0.6) +
    geom_line(data = avg, aes(group = 1), colour = MEAN_COL,
              linewidth = 0.9) +
    geom_point(data = avg, colour = MEAN_COL, size = 1.3) +
    geom_point(data = dg, shape = 21, size = 3, stroke = 1,
               fill = NA, colour = REF_COL) +
    facet_grid(. ~ ct) +
    labs(x = if (xaxis) "c Used for Fitting" else NULL, y = ylab) +
    theme_cowplot(font_size = 11) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "grey80", fill = NA),
          axis.title.y = element_text(size = 10.5),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          strip.text = if (strips) element_text(face = "bold", size = 10.5)
                       else element_blank())
  if (logy) g <- g + scale_y_log10()
  g
}

truth_rc <- aggregate(ex$rankcor_F_abs_true, by = list(ct = ex$ct), FUN = mean)
names(truth_rc)[2] <- "value"

pA <- row_panel("rel_rmse", "Rel. RMSE of Rates", ref1 = TRUE, logy = TRUE,
                strips = TRUE)
pB <- row_panel("rankcor_F_abs", "Mean Abs. Factor Correlation",
                truth_df = truth_rc, xaxis = TRUE)

g <- plot_grid(pA, pB, ncol = 1, align = "v", axis = "lr",
               labels = c("A", "B"), label_size = 13,
               rel_heights = c(1.06, 1.04))

ggsave("../images/sim2_rel_rmse.png", g, device = "png",
       width = 8.5, height = 6.2, dpi = 300, bg = "white")
message("Wrote ../images/sim2_rel_rmse.png")
