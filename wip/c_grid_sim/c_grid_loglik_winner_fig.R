# Which c_fit attains the highest log-likelihood, as a function of the c that
# generated the data? Counts over 10 seeds, as an 8 x 8 grid with the diagonal
# (c_fit == c_true) outlined.
#
# Reads c_grid_sim_metrics.rds (written by experiments/analyze_c_grid_sim.R).
#
# Comparing log-likelihoods across c_fit is legitimate: every fit scores the
# same y under the same dominating measure with constants included, so no
# Jacobian enters. Each cell uses the better of its two inits by
# log-likelihood, though that is immaterial here -- the two agree to within
# 0.008 on all 640 cells.
#
# c is a discrete axis rather than a numeric one because c_fit = Inf is a real
# level of the grid (the c -> Inf limit, plain Poisson NMF) and has no place on
# a log axis.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))

exp_dir <- "../experiments"
out_png <- "../images/c_grid_loglik_winner.png"

d <- readRDS(file.path(exp_dir, "c_grid_sim_metrics.rds"))

## best init per (seed, c_true, c_fit)
key <- paste(d$seed, d$c_true, d$c_fit)
b   <- d[order(key, -d$loglik), ]
b   <- b[!duplicated(paste(b$seed, b$c_true, b$c_fit)), ]

CLEV <- sort(unique(b$c_fit))
CLAB <- ifelse(is.infinite(CLEV), "∞",
               format(CLEV, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
fc <- function(x) factor(x, levels = CLEV, labels = CLAB)

win <- do.call(rbind, lapply(split(b, paste(b$seed, b$c_true)), function(x)
  x[which.max(x$loglik), c("seed", "c_true", "c_fit")]))

cnt <- as.data.frame(table(ct = fc(win$c_true), cf = fc(win$c_fit)),
                     responseName = "n")
cnt$lab <- ifelse(cnt$n == 0, "", cnt$n)
diagd <- data.frame(ct = CLAB, cf = CLAB)

p <- ggplot(cnt, aes(cf, ct)) +
  geom_tile(aes(fill = n), colour = "white", linewidth = 0.6) +
  geom_tile(data = diagd, fill = NA, colour = "#C2410C", linewidth = 0.9) +
  geom_text(aes(label = lab, colour = n > 5), size = 4, show.legend = FALSE) +
  scale_fill_gradient(low = "#F2F6F7", high = "#0E5C6B", name = "seeds\nwon",
                      limits = c(0, 10), breaks = c(0, 5, 10)) +
  scale_colour_manual(values = c("FALSE" = "grey15", "TRUE" = "white")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "c used for fitting", y = "c used to simulate",
       title = "Which c achieves the highest log-likelihood",
       subtitle = "counts over 10 seeds; outlined cells are c_fit = c_true") +
  coord_fixed() +
  theme(plot.title = element_text(size = 13, face = "plain"),
        plot.subtitle = element_text(size = 10, colour = "grey35"),
        panel.border = element_rect(colour = "grey70", fill = NA))

ggsave(out_png, p, width = 7, height = 6, dpi = 150, bg = "white")
message("Wrote ", normalizePath(out_png, mustWork = FALSE))
