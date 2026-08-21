# Appendix figure: the peak-to-mean ratio of the fitted representation
# pm(B), B = L F', for every exact fit of Simulation 2, against the value
# PREDICTED from the true rates alone: pm(log(1 + Lambda/c)) at finite
# fitting c, and pm(Lambda) at c = Inf. Supports the main-text claim that
# the fitted pm(B) closely tracks the prediction from Lambda even when
# Lambda was generated at a different c_true.
#
# One row: fitted pm(B) (thin lines: replicates; thick: mean) with the
# mean prediction from Lambda dashed, and -- as a contrast -- the DATA
# budget pm(log(1 + Y/c)) (pm(Y) at c = Inf) dot-dashed: the version
# computable without knowing the truth, which the zeros inflate at small
# c.
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
    d <- sim_dataset(s, ct)
    data.frame(seed = s, c_true = ct, cc = CC,
               pm_pred = sapply(CC, pm_of, M = d$lambda),
               pm_data = sapply(CC, pm_of, M = d$Y))
  }))
}))
message("predictions computed for ", nrow(pred), " (seed, c_true, c) triples")

key <- function(d) paste(d$seed, format(d$c_true), format(d$cc))
ex$pm_pred <- pred$pm_pred[match(key(ex), key(pred))]
ex$pm_data <- pred$pm_data[match(key(ex), key(pred))]
stopifnot(!anyNA(ex$pm_pred), !anyNA(ex$pm_data))

SEED_COL <- "#9BB7BD"; MEAN_COL <- "#0E5C6B"
PRED_COL <- "#C2410C"; DATA_COL <- "#B45309"

avg <- aggregate(cbind(pm_B, pm_pred, pm_data) ~ cf + ct, ex, mean)
g <- ggplot(ex, aes(cf, pm_B)) +
  geom_line(aes(group = seed), colour = SEED_COL,
            linewidth = 0.3, alpha = 0.6) +
  geom_line(data = avg, aes(group = 1), colour = MEAN_COL, linewidth = 0.9) +
  geom_point(data = avg, colour = MEAN_COL, size = 1.3) +
  geom_line(data = avg, aes(y = pm_pred, group = 1), colour = PRED_COL,
            linetype = "dashed", linewidth = 0.7) +
  geom_line(data = avg, aes(y = pm_data, group = 1), colour = DATA_COL,
            linetype = "dotdash", linewidth = 0.7) +
  scale_y_log10(labels = function(x)
    format(x, scientific = FALSE, big.mark = ",", trim = TRUE)) +
  facet_grid(. ~ ct) +
  labs(x = "c Used for Fitting", y = "pm(B) of Fitted Model") +
  theme_cowplot(font_size = 11) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "grey80", fill = NA),
        axis.title.y = element_text(size = 10.5),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 10.5))

ggsave("../images/sim2_pm_budget.png", g, device = "png",
       width = 8.5, height = 3.4, dpi = 300, bg = "white")
message("Wrote ../images/sim2_pm_budget.png")
