# Which c_fit maximizes the log-likelihood, as a function of the c that
# generated the data?
#
#   A  winner distribution: how often each c_fit wins, over 10 seeds, for each
#      c_true (counts in an 8 x 8 grid, diagonal outlined)
#   B  mean suboptimality vs c_fit, one curve per c_true
#
# Reads c_grid_sim_metrics.rds (written by experiments/analyze_c_grid_sim.R).
#
# Suboptimality for a fit = (best log-likelihood any c_fit achieved on THAT
# seed's dataset) - (this fit's log-likelihood). It is 0 for the winner and
# positive elsewhere, and it is a fair comparison because every c_fit is scored
# on the same y under the same dominating measure -- these are all Poisson
# models for the same counts, constants included, so no Jacobian enters.
#
# Each cell was fit from a rank-1 and a random start; we take the better of the
# two by log-likelihood. In this simulation that choice is immaterial -- the two
# agree to within 0.008 everywhere -- but it is the right summary for a
# nonconvex objective.
#
# Both axes carry c on a discrete scale rather than a numeric one, because
# c_fit = Inf is a real level of the grid (the c -> Inf limit, plain Poisson
# NMF) and cannot be placed on a log axis. y in panel B is log1p: suboptimality
# is exactly 0 at the winner, which log10 could not show, and differences below
# ~1 log-likelihood unit are meaningless anyway.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))

exp_dir <- "../experiments"
out_png <- "../images/c_grid_loglik.png"

d <- readRDS(file.path(exp_dir, "c_grid_sim_metrics.rds"))

## ---- best init per (seed, c_true, c_fit) -----------------------------------
key <- paste(d$seed, d$c_true, d$c_fit)
b   <- d[order(key, -d$loglik), ]
b   <- b[!duplicated(paste(b$seed, b$c_true, b$c_fit)), ]

## suboptimality within each (seed, c_true)
b$subopt <- ave(b$loglik, paste(b$seed, b$c_true), FUN = max) - b$loglik

## ---- shared discrete c scale -----------------------------------------------
CLEV <- sort(unique(b$c_fit))                       # 1e-3 ... 1e3, Inf
CLAB <- ifelse(is.infinite(CLEV), "∞",
               format(CLEV, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
fc <- function(x) factor(x, levels = CLEV, labels = CLAB)

b$cf <- fc(b$c_fit)
b$ct <- fc(b$c_true)

## ---- A: winner distribution ------------------------------------------------
win <- do.call(rbind, lapply(split(b, paste(b$seed, b$c_true)), function(x)
  x[which.max(x$loglik), c("seed", "c_true", "c_fit")]))

cnt <- as.data.frame(table(ct = fc(win$c_true), cf = fc(win$c_fit)),
                     responseName = "n")
cnt$lab <- ifelse(cnt$n == 0, "", cnt$n)
diagd <- data.frame(ct = CLAB, cf = CLAB)          # where c_fit == c_true

pA <- ggplot(cnt, aes(cf, ct)) +
  geom_tile(aes(fill = n), colour = "white", linewidth = 0.6) +
  geom_tile(data = diagd, fill = NA, colour = "#C2410C", linewidth = 0.9) +
  geom_text(aes(label = lab, colour = n > 5), size = 3.6, show.legend = FALSE) +
  scale_fill_gradient(low = "#F2F6F7", high = "#0E5C6B", name = "seeds\nwon",
                      limits = c(0, 10), breaks = c(0, 5, 10)) +
  scale_colour_manual(values = c("FALSE" = "grey15", "TRUE" = "white")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "c used for fitting", y = "c used to simulate",
       title = "A  Which c wins, out of 10 seeds") +
  coord_fixed() +
  theme(plot.title = element_text(size = 12, face = "plain"),
        panel.border = element_rect(colour = "grey70", fill = NA))

## ---- B: mean suboptimality -------------------------------------------------
## sequential ramp: c_true is ordered, so a single-hue light-to-dark ramp is the
## right encoding rather than categorical hues. Direct labels at the right edge
## carry the identity, so the ramp only has to convey order.
agg <- aggregate(subopt ~ ct + cf, data = b, FUN = mean)
## every step clears 3:1 against white (min 3.37) with monotone decreasing
## lightness -- a paler start looked better but the first two steps were
## invisible as lines. Direct labels carry identity, so the ramp only has to
## convey order.
RAMP <- colorRampPalette(c("#4E96A5", "#06333C"))(length(CLAB))
names(RAMP) <- CLAB

ends <- agg[agg$cf == tail(CLAB, 1), ]

pB <- ggplot(agg, aes(cf, subopt, colour = ct, group = ct)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  ## mark where the fitted c equals the simulating c
  geom_point(data = agg[as.character(agg$ct) == as.character(agg$cf), ],
             shape = 21, size = 3.4, stroke = 1, fill = NA, colour = "#C2410C") +
  ggrepel::geom_text_repel(data = ends, aes(label = ct), nudge_x = 0.6,
                           direction = "y", hjust = 0, size = 3.2, min.segment.length = 0,
                           segment.colour = "grey75", segment.size = 0.3,
                           box.padding = 0.12, show.legend = FALSE) +
  scale_colour_manual(values = RAMP, name = "c simulated") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000),
                     labels = c("0", "1", "10", "100", "1000")) +
  expand_limits(x = length(CLAB) + 1.6) +
  labs(x = "c used for fitting",
       y = "mean log-likelihood below the best c",
       caption = "circled: c_fit = c_true",
       title = "B  Cost of fitting at the wrong c") +
  theme(plot.title = element_text(size = 12, face = "plain"),
        plot.caption = element_text(size = 9, colour = "grey35", hjust = 1),
        legend.position = "none")

final <- plot_grid(pA, pB, nrow = 1, rel_widths = c(1, 1.15), align = "h", axis = "tb")
ggsave(out_png, final, width = 12, height = 5.2, dpi = 150, bg = "white")
message("Wrote ", normalizePath(out_png, mustWork = FALSE))
