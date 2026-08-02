# Log-likelihood across fitted c, one panel per simulating c. Each panel shows
# all 10 seeds individually plus their mean, so the seed-to-seed spread is
# visible rather than averaged away.
#
# Reads c_grid_sim_metrics.rds (written by experiments/analyze_c_grid_sim.R).
#
# y is the log-likelihood measured DOWN FROM that seed's best c, i.e.
#   max_over_c_fit(loglik for this seed) - loglik.
# Raw log-likelihood cannot be overlaid across seeds: each seed is a different
# dataset, so the curves would be offset by dataset-level constants of a few
# hundred units and the shape -- the thing being compared -- would be
# unreadable. Subtracting each seed's own maximum removes exactly that constant
# and leaves 0 at the winner.
#
# The comparison across c_fit within a seed is fair: all eight are Poisson
# models for the same y under the same dominating measure with constants
# included, so no Jacobian enters.
#
# y is log1p, because suboptimality is exactly 0 at the winner (which log10
# cannot show) and differences below ~1 log-likelihood unit are meaningless.
# x is discrete because c_fit = Inf is a real level of the grid.
#
# Panels share a y scale on purpose -- the point is that the penalty for being
# in the wrong REGIME (~1000+) dwarfs the penalty within a regime (~1), and a
# free y would hide that by rescaling every panel to its own range.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 11))

exp_dir <- "../experiments"
out_png <- "../images/c_grid_loglik_by_ctrue.png"

d <- readRDS(file.path(exp_dir, "c_grid_sim_metrics.rds"))

## best init per (seed, c_true, c_fit); the two agree to within 0.008 here
key <- paste(d$seed, d$c_true, d$c_fit)
b   <- d[order(key, -d$loglik), ]
b   <- b[!duplicated(paste(b$seed, b$c_true, b$c_fit)), ]
b$subopt <- ave(b$loglik, paste(b$seed, b$c_true), FUN = max) - b$loglik

CLEV <- sort(unique(b$c_fit))
CLAB <- ifelse(is.infinite(CLEV), "∞",
               format(CLEV, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
fc <- function(x) factor(x, levels = CLEV, labels = CLAB)
b$cf <- fc(b$c_fit)
b$ct <- factor(paste("simulated at c =", fc(b$c_true)),
               levels = paste("simulated at c =", CLAB))

## per-panel mean over seeds, and the position of the matching c_fit
avg  <- aggregate(subopt ~ ct + cf, data = b, FUN = mean)
diag <- unique(b[as.character(b$cf) == as.character(fc(b$c_true)),
                 c("ct", "cf")])
diag <- merge(diag, avg, by = c("ct", "cf"))

p <- ggplot(b, aes(cf, subopt)) +
  ## individual seeds, kept visible so the spread is not hidden by the mean
  geom_line(aes(group = seed), colour = "#9BB7BD", linewidth = 0.35, alpha = 0.8) +
  geom_line(data = avg, aes(group = 1), colour = "#0E5C6B", linewidth = 1) +
  geom_point(data = avg, colour = "#0E5C6B", size = 1.6) +
  ## where the fitted c equals the simulating c
  geom_point(data = diag, shape = 21, size = 3.6, stroke = 1.1,
             fill = NA, colour = "#C2410C") +
  facet_wrap(~ ct, nrow = 2) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000),
                     labels = c("0", "1", "10", "100", "1000")) +
  labs(x = "c used for fitting",
       y = "log-likelihood below that seed's best c",
       title = "Cost of fitting at the wrong c",
       subtitle = paste("thin lines: 10 individual seeds;  thick line: mean;",
                        "circled: c_fit = c_true")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 13, face = "plain"),
        plot.subtitle = element_text(size = 10, colour = "grey35"),
        panel.border = element_rect(colour = "grey80", fill = NA),
        axis.text.x = element_text(size = 8))

ggsave(out_png, p, width = 12, height = 6.4, dpi = 150, bg = "white")
message("Wrote ", normalizePath(out_png, mustWork = FALSE))
