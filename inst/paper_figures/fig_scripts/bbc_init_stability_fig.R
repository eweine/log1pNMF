# BBC initialization-stability experiment (reviewer request): exact
# log-likelihood of the converged K = 10 BBC fits across c in
# {0.001, 1, Inf}, from 30 random initializations plus the paper's rank-1
# initialization.
#
# One row: a panel of the log-likelihoods centered at each c's random-init
# mean (individual seeds, +/- 1 sd bars, and the rank-1 fit), and a table
# of the raw mean and sd of the random-init log-likelihoods per c. The
# centering is what makes the initialization spread visible: the between-c
# differences (~60k log-lik units) are ~40x the within-c sd.
#
# Reads res_list$bbc_init_stability from ../data/experiment_results.Rdata
# (the summary rows of the per-task outputs of run_bbc_init.R, collected by
# utils/collect_data.R).

library(ggplot2)
library(cowplot)

load("../data/experiment_results.Rdata")
rows <- res_list$bbc_init_stability
stopifnot(nrow(rows) == 93, all(rows$converged))

fmt_c <- function(x) ifelse(is.infinite(x), "∞",
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
rows$cf <- factor(fmt_c(rows$cc), levels = fmt_c(sort(unique(rows$cc))))

rnd <- subset(rows, init == "random")
r1  <- subset(rows, init == "rank1")
agg <- do.call(rbind, lapply(split(rnd, rnd$cf), function(d)
  data.frame(cf = d$cf[1], mean = mean(d$loglik_exact),
             sd = sd(d$loglik_exact))))
r1$center  <- r1$loglik_exact  - agg$mean[match(r1$cf, agg$cf)]
rnd$center <- rnd$loglik_exact - agg$mean[match(rnd$cf, agg$cf)]

MEAN_COL <- "#0E5C6B"; SEED_COL <- "#9BB7BD"; R1_COL <- "#C2410C"

g <- ggplot(rnd, aes(cf, center)) +
  geom_jitter(width = 0.08, height = 0, colour = SEED_COL,
              size = 1.4, alpha = 0.8) +
  geom_errorbar(data = agg, aes(cf, y = 0, ymin = -sd, ymax = sd),
                width = 0.12, colour = MEAN_COL, linewidth = 0.6) +
  geom_point(data = agg, aes(cf, 0, colour = "Random inits"), size = 2.4) +
  geom_point(data = r1, aes(cf, center, colour = "Rank-1 init"),
             shape = 18, size = 3) +
  scale_colour_manual(values = c("Random inits" = MEAN_COL,
                                 "Rank-1 init" = R1_COL), name = NULL,
                      guide = guide_legend(
                        override.aes = list(shape = c(16, 18), size = 3))) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "c Used for Fitting",
       y = "Log-Likelihood − Random-Init Mean",
       title = "BBC K = 10 Replicates") +
  theme_cowplot(font_size = 12) +
  theme(panel.border = element_rect(colour = "grey80", fill = NA),
        plot.title = element_text(size = 13, hjust = 0.5),
        legend.position = "right")

ggsave("../images/bbc_init_stability.png", g, width = 7, height = 4.2,
       dpi = 300, bg = "white")
message("Wrote ../images/bbc_init_stability.png")

## LaTeX table companion (pasted into the manuscript's subfigure): raw mean
## and sd of the random-init exact log-likelihoods per c.
cat("\nLaTeX table rows:\n")
for (i in seq_len(nrow(agg)))
  cat(sprintf("%s & $%s$ & $%s$ \\\\\n",
              ifelse(agg$cf[i] == "∞", "$\\infty$", paste0("$", agg$cf[i], "$")),
              format(round(agg$mean[i]), big.mark = ","),
              format(round(agg$sd[i]), big.mark = ",")))
