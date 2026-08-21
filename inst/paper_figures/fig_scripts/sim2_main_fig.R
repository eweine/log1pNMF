# Main figure for Simulation 2: exact fits across the fitting-c grid
# (including c = Inf, the topic model) on data generated at
# c_true in {1e-3, 1, Inf}. Four rows x three columns (one per c_true):
#
#   A  mean Hoyer sparsity of the loadings L
#   B  mean Hoyer sparsity of the factors F
#   C  per-entry likelihood ratio vs the fit at c_fit = c_true
#      (exp of the mean per-entry exact log-likelihood difference)
#   D  RMSE of the fitted rates against the true rates, relative to the
#      fit at c_fit = c_true (log scale)
#
# Thin lines: individual seeds; thick line: mean over seeds; dashed: the
# truth's own value (A-B) or the reference ratio 1 (C-D); open circle:
# c_fit = c_true. Reads approx_sim_metrics.rds (analyze_approx_sim.R);
# rerunning after more seeds land just adds thin lines.

library(ggplot2)
library(cowplot)

res <- readRDS("../experiments/approx_sim_metrics.rds")
ex  <- subset(res, method == "exact")

NP <- unique(res$n * res$p)
stopifnot(length(NP) == 1)

fmt_c <- function(x) ifelse(is.infinite(x), "∞",
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
CC  <- sort(unique(ex$cc))                      # Inf sorts last
ex$cf <- factor(fmt_c(ex$cc), levels = fmt_c(CC))
CTL <- paste0("Data Simulated with c = ", fmt_c(sort(unique(ex$c_true))))
ex$ct <- factor(paste0("Data Simulated with c = ", fmt_c(ex$c_true)),
                levels = CTL)
ex$is_ref <- (is.finite(ex$c_true) & ex$cc == ex$c_true) |
             (is.infinite(ex$c_true) & is.infinite(ex$cc))

## references: the same seed's fit at c_fit = c_true
ref <- ex[ex$is_ref, c("seed", "ct", "loglik_exact", "rmse")]
names(ref)[3:4] <- c("ll_ref", "rmse_ref")
ex <- merge(ex, ref, by = c("seed", "ct"))
ex$lr_pe    <- exp((ex$loglik_exact - ex$ll_ref) / NP)
ex$rel_rmse <- ex$rmse / ex$rmse_ref

## ---- one row of the figure --------------------------------------------------
SEED_COL <- "#9BB7BD"; MEAN_COL <- "#0E5C6B"; REF_COL <- "#C2410C"

row_panel <- function(d, ylab, truth_df = NULL, ref1 = FALSE, logy = FALSE,
                      strips = FALSE, xaxis = FALSE) {
  avg <- aggregate(value ~ cf + ct, d, mean)
  dg  <- merge(unique(d[d$is_ref, c("cf", "ct")]), avg)
  p <- ggplot(d, aes(cf, value))
  if (ref1)
    p <- p + geom_hline(yintercept = 1, linetype = "dashed",
                        colour = REF_COL, linewidth = 0.4)
  if (!is.null(truth_df))
    p <- p + geom_hline(data = truth_df, aes(yintercept = value),
                        linetype = "dashed", colour = REF_COL,
                        linewidth = 0.4)
  p <- p +
    geom_line(aes(group = seed), colour = SEED_COL,
              linewidth = 0.3, alpha = 0.6) +
    geom_line(data = avg, aes(group = 1), colour = MEAN_COL,
              linewidth = 0.9) +
    geom_point(data = avg, colour = MEAN_COL, size = 1.3) +
    geom_point(data = dg, shape = 21, size = 3, stroke = 1,
               fill = NA, colour = REF_COL) +
    facet_grid(. ~ ct) +
    labs(x = if (xaxis) "c used for fitting" else NULL, y = ylab) +
    theme_cowplot(font_size = 11) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "grey80", fill = NA),
          axis.title.y = element_text(size = 10.5),
          strip.text = if (strips) element_text(face = "bold", size = 10.5)
                       else element_blank())
  if (logy) p <- p + scale_y_log10()
  p <- p + theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
  p
}

long <- function(v) { d <- ex; d$value <- ex[[v]]; d }
truth_of <- function(v) {
  td <- aggregate(ex[[v]], by = list(ct = ex$ct), FUN = mean)
  names(td)[2] <- "value"; td
}

pA <- row_panel(long("hoyer_L"), "Mean Loading Sparsity",
                truth_df = truth_of("hoyer_L_true"), strips = TRUE)
pB <- row_panel(long("hoyer_F"), "Mean Factor Sparsity",
                truth_df = truth_of("hoyer_F_true"))
pC <- row_panel(long("lr_pe"), "Likelihood Ratio", ref1 = TRUE)
pD <- row_panel(long("rel_rmse"), "Relative RMSE of Rates", ref1 = TRUE,
                logy = TRUE, xaxis = TRUE)

g <- plot_grid(pA, pB, pC, pD, ncol = 1, align = "v", axis = "lr",
               labels = c("A", "B", "C", "D"), label_size = 13,
               rel_heights = c(1.1, 1, 1, 1.04))

ggsave("../images/sim2_main.png", g, device = "png",
       width = 8, height = 10, dpi = 300, bg = "white")
message("Wrote ../images/sim2_main.png")
