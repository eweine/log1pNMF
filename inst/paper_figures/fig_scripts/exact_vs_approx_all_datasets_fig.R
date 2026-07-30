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
# that dataset) - (current exact log-likelihood), on a log axis. The REFERENCE is
# still a single per-dataset value, so the four curves of a dataset are all
# measured against the same target; only the visible window differs per panel.
#
# Every panel gets its own x AND y scale (facet scales = "free"). The panels
# within a row span very different ranges -- pancreas rank-1 reaches ~1e0 while
# pancreas random plateaus at ~1e5 -- and forcing a shared y flattened both into
# a decade-wide sliver. The cost is that you can no longer read relative height
# across the two panels of a row by eye; compare the axis labels, not the curves.
#
# Time uses a LOG1P axis for pancreas and BBC (xtrans below). Both essentially
# converge in their first seconds and then flatline, so on a linear axis all the
# structure was crushed against x = 0. A plain log axis cannot be used: iterate 0
# is the INITIALIZATION, recorded at t = 0 before any update (see the pre-loop
# loglik_history.push_back in pois_reg_log1p.cpp, and the matching t = 0 seed in
# exact_vs_approx_helpers.R). That point matters -- under the random scheme both
# methods start from one shared, identical initialization, which is the anchor
# the whole comparison rests on -- and log10 would silently drop it, making the
# curves appear to start at different qualities. log1p keeps t = 0 at 0.
#
# Those two rows are therefore in SECONDS, not hours/minutes: log1p(x) is ~x
# below 1 and ~log(x) above, so the knee sits at x = 1 in whatever unit is used.
# In hours the entire interesting range (0.4 s to a few minutes) falls below 1
# and stays compressed; in seconds it spreads out properly.
#
# MOCA stays linear in hours: its curves improve steadily across the whole 10 h
# budget rather than converging early, so a log axis would misrepresent them.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

exp_dir <- "../experiments"
out_png <- "../images/exact_vs_approx_all_datasets.png"

## dataset tag, display title, wall-time unit (divisor + label), and whether the
## time axis is log1p-transformed
DATASETS <- list(
  list(key = "moca",     tag = "exact_vs_approx_moca_K25_c1",
       title = "MOCA — 1.33M cells × 26,182 genes (K = 25)",
       tdiv = 3600, tlab = "Wall-clock time (hours)",   xlog1p = FALSE),
  list(key = "pancreas", tag = "exact_vs_approx_pancreas_K13_c1",
       title = "Pancreas cytokine (K = 13)",
       tdiv = 1,    tlab = "Wall-clock time (seconds)", xlog1p = TRUE),
  list(key = "bbc",      tag = "exact_vs_approx_bbc_K10_c1",
       title = "BBC news — documents × terms (K = 10)",
       tdiv = 1,    tlab = "Wall-clock time (seconds)", xlog1p = TRUE)
)

## plain decimal labels ("0", "1", "10", "1000") rather than scientific notation
dec_lab <- function(x) ifelse(is.na(x), "",
                              formatC(x, format = "fg", digits = 3, drop0trailing = TRUE))

## 0 plus the powers of ten spanned by the panel -- log1p is linear below 1, so
## a decade grid starting at 0 reads naturally
log1p_breaks <- function(lims) {
  hi <- max(lims, na.rm = TRUE)
  b  <- c(0, 10^(0:ceiling(log10(max(hi, 1)))))
  b[b <= hi]
}

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
  ## log1p keeps t = 0, so the initialization iterate stays in the plot
  d$hours  <- d$time / ds$tdiv
  d$scheme <- factor(SCHEME_LABS[d$scheme], levels = SCHEME_LABS)
  d$method <- factor(ifelse(d$method == "approx", "Approximate", "Exact"),
                     levels = c("Approximate", "Exact"))
  g <- ggplot(d, aes(hours, excess, color = method)) +
    geom_line(linewidth = 0.8) +
    ## both scales free per panel -- see the header note on the tradeoff
    facet_wrap(~ scheme, nrow = 1, scales = "free") +
    scale_color_manual(values = MCOL, name = "Optimization") +
    scale_y_log10() +
    labs(x = ds$tlab, y = NULL, title = ds$title) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(size = 12, face = "plain"),
          legend.position = "none")
  if (ds$xlog1p)
    g <- g + scale_x_continuous(trans = "log1p", breaks = log1p_breaks,
                                labels = dec_lab)
  g
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
