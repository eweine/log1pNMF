# Appendix figure: the peak-to-mean ratio of the fitted representation
# pm(B), B = L F', for every exact fit of Simulation 2, against the value
# PREDICTED from the true rates alone: pm(log(1 + Lambda/c)) at finite
# fitting c, and pm(Lambda) at c = Inf. Supports the main-text claim that
# the fitted pm(B) closely tracks the prediction from Lambda even when
# Lambda was generated at a different c_true.
#
# Top row: fitted pm(B) (thin lines: replicates; thick: mean) with the
# mean prediction dashed. Bottom row: the ratio fitted / predicted, per
# replicate, with the mean; dashed line at 1.
#
# Reads approx_sim_metrics.rds and regenerates the true rates from the
# deterministic simulation generator (approx_sim_common.R), so it needs no
# additional downloaded outputs and picks up extra seeds automatically.

library(ggplot2)
library(cowplot)

res <- readRDS("../experiments/approx_sim_metrics.rds")
source("../experiments/approx_sim_common.R", chdir = TRUE)
ex  <- subset(res, method == "exact")

fmt_c <- function(x) ifelse(is.infinite(x), "∞",
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
CC <- sort(unique(ex$cc))
ex$cf <- factor(fmt_c(ex$cc), levels = fmt_c(CC))
CTL <- paste0("Data Simulated with c = ", fmt_c(sort(unique(ex$c_true))))
ex$ct <- factor(paste0("Data Simulated with c = ", fmt_c(ex$c_true)),
                levels = CTL)

## predicted budget from the true rates, per (c_true, seed, fitting c)
pm_of <- function(M, cc) {
  g <- if (is.finite(cc)) log1p(M / cc) else M
  max(g) / mean(g)
}
pred <- do.call(rbind, lapply(sort(unique(ex$c_true)), function(ct) {
  seeds <- sort(unique(ex$seed[
    (is.finite(ex$c_true) & ex$c_true == ct) |
    (is.infinite(ex$c_true) & is.infinite(ct))]))
  do.call(rbind, lapply(seeds, function(s) {
    lam <- sim_dataset(s, ct)$lambda
    data.frame(seed = s, c_true = ct, cc = CC,
               pm_pred = sapply(CC, pm_of, M = lam))
  }))
}))
message("predictions computed for ", nrow(pred), " (seed, c_true, c) triples")

key <- function(d) paste(d$seed, format(d$c_true), format(d$cc))
ex$pm_pred <- pred$pm_pred[match(key(ex), key(pred))]
stopifnot(!anyNA(ex$pm_pred))
ex$ratio <- ex$pm_B / ex$pm_pred

SEED_COL <- "#9BB7BD"; MEAN_COL <- "#0E5C6B"; REF_COL <- "#C2410C"

theme_row <- function(strips) theme(
  strip.background = element_blank(),
  panel.border = element_rect(colour = "grey80", fill = NA),
  axis.title.y = element_text(size = 10.5),
  axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
  strip.text = if (strips) element_text(face = "bold", size = 10.5)
               else element_blank())

avg <- aggregate(cbind(pm_B, pm_pred) ~ cf + ct, ex, mean)
pA <- ggplot(ex, aes(cf, pm_B)) +
  geom_line(aes(group = seed), colour = SEED_COL,
            linewidth = 0.3, alpha = 0.6) +
  geom_line(data = avg, aes(group = 1), colour = MEAN_COL, linewidth = 0.9) +
  geom_point(data = avg, colour = MEAN_COL, size = 1.3) +
  geom_line(data = avg, aes(y = pm_pred, group = 1), colour = REF_COL,
            linetype = "dashed", linewidth = 0.7) +
  scale_y_log10(labels = function(x)
    format(x, scientific = FALSE, big.mark = ",", trim = TRUE)) +
  facet_grid(. ~ ct) +
  labs(x = NULL, y = "pm(B) of Fitted Model") +
  theme_cowplot(font_size = 11) + theme_row(strips = TRUE)

ravg <- aggregate(ratio ~ cf + ct, ex, mean)
pB <- ggplot(ex, aes(cf, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL,
             linewidth = 0.4) +
  geom_line(aes(group = seed), colour = SEED_COL,
            linewidth = 0.3, alpha = 0.6) +
  geom_line(data = ravg, aes(group = 1), colour = MEAN_COL, linewidth = 0.9) +
  geom_point(data = ravg, colour = MEAN_COL, size = 1.3) +
  scale_y_log10(breaks = c(0.5, 0.707, 1, 1.414),
                labels = c("0.5", "0.71", "1", "1.41")) +
  facet_grid(. ~ ct) +
  labs(x = "c Used for Fitting", y = "Fitted / Predicted") +
  theme_cowplot(font_size = 11) + theme_row(strips = FALSE)

g <- plot_grid(pA, pB, ncol = 1, align = "v", axis = "lr",
               labels = c("A", "B"), label_size = 13,
               rel_heights = c(1.16, 1))

ggsave("../images/sim2_pm_budget.png", g, device = "png",
       width = 8, height = 5.4, dpi = 300, bg = "white")
message("Wrote ../images/sim2_pm_budget.png")
