# Distance-from-best-log-likelihood vs. wall-clock time, approximate vs. exact
# optimization, as a 3 x 2 grid: one row per dataset (MOCA, pancreas, BBC), one
# column per initialization scheme (random, rank-1).
#
# Reads the per-(scheme, method) trace tables written by the experiment jobs
# (exact_vs_approx_<dataset>_..._traces.rds), each a tidy data frame with columns
# time (seconds), exact_loglik, scheme, method, iter. For the rank-1 scheme the
# method-appropriate warm-up is already chained in at the front of each curve.
#
# y = distance from the best log-likelihood = (best exact log-likelihood found for
# that dataset) - (current exact log-likelihood), on a log axis. A single per-
# dataset reference keeps the two schemes of a dataset directly comparable. The
# datasets live on very different time scales (MOCA/pancreas ~10 h, BBC seconds-
# minutes), so each row gets its own wall-time axis.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

exp_dir <- "../experiments"
out_png <- "../images/exact_vs_approx_all_datasets.png"

## dataset tag, display title, and wall-time unit (divisor + label)
DATASETS <- list(
  list(key = "moca",     tag = "exact_vs_approx_moca_K25_c1",
       title = "MOCA — 1.33M cells × 26,182 genes (K = 25)",
       tdiv = 3600, tlab = "Wall-clock time (hours)"),
  list(key = "pancreas", tag = "exact_vs_approx_pancreas_K13_c1",
       title = "Pancreas cytokine (K = 13)",
       tdiv = 3600, tlab = "Wall-clock time (hours)"),
  list(key = "bbc",      tag = "exact_vs_approx_bbc_K10_c1",
       title = "BBC news — documents × terms (K = 10)",
       tdiv = 60,   tlab = "Wall-clock time (minutes)")
)

SCHEME_LABS <- c(random = "Random init", rank1 = "Rank-1 init")
## approx vs exact -- validated colorblind-safe pair (matches the other timing figs)
MCOL <- c("Approximate" = "#E69F00", "Exact" = "#0072B2")

load_dataset <- function(tag) {
  combos <- expand.grid(scheme = c("random", "rank1"), method = c("approx", "exact"),
                        stringsAsFactors = FALSE)
  d <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
    f <- file.path(exp_dir, sprintf("%s_%s_%s_traces.rds", tag, combos$scheme[i], combos$method[i]))
    if (!file.exists(f)) { message("  missing: ", basename(f)); return(NULL) }
    readRDS(f)
  }))
  d[is.finite(d$exact_loglik), ]
}

## one dataset -> a two-panel (random | rank-1) ggplot, no legend
panel_for <- function(ds) {
  d <- load_dataset(ds$tag)
  best     <- max(d$exact_loglik)                 # single reference for this dataset
  d$excess <- best - d$exact_loglik
  d        <- d[d$excess > 0, ]                    # drop the point that attains best
  d$hours  <- d$time / ds$tdiv
  d$scheme <- factor(SCHEME_LABS[d$scheme], levels = SCHEME_LABS)
  d$method <- factor(ifelse(d$method == "approx", "Approximate", "Exact"),
                     levels = c("Approximate", "Exact"))
  ggplot(d, aes(hours, excess, color = method)) +
    geom_line(linewidth = 0.8) +
    ## x free per panel (datasets/schemes converge on very different time scales);
    ## y shared within the dataset so the two schemes stay directly comparable
    facet_wrap(~ scheme, nrow = 1, scales = "free_x") +
    scale_color_manual(values = MCOL, name = "Optimization") +
    scale_y_log10() +
    labs(x = ds$tlab, y = NULL, title = ds$title) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(size = 12, face = "plain"),
          legend.position = "none")
}

panels <- lapply(DATASETS, panel_for)

## shared legend (pull from a throwaway plot that has one)
legend <- get_legend(panels[[1]] + theme(legend.position = "top"))

grid  <- plot_grid(plotlist = panels, ncol = 1, align = "v", axis = "lr")
## a single rotated y-axis label shared by all rows
ylab  <- ggdraw() + draw_label("Distance from best log-likelihood",
                               angle = 90, size = 13)
body  <- plot_grid(ylab, grid, ncol = 2, rel_widths = c(0.04, 1))
final <- plot_grid(legend, body, ncol = 1, rel_heights = c(0.05, 1))

ggsave(out_png, final, width = 9, height = 11, dpi = 150, bg = "white")
message("Wrote ", normalizePath(out_png, mustWork = FALSE))
