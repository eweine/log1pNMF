# Distance-from-best-log-likelihood vs. wall-clock time, approximate vs. exact
# optimization, as one row of three panels: one per dataset (MOCA, pancreas,
# BBC), RANK-1 INITIALIZATION ONLY (the paper's default; the random-init
# comparison lives in the git history of this script).
#
# Reads the per-(scheme, method) trace tables written by the experiment jobs
# (exact_vs_approx_<dataset>_..._traces.rds), each a tidy data frame with columns
# time (seconds), exact_loglik, scheme, method, iter. The method-appropriate
# rank-1 warm-up is already chained in at the front of each curve.
#
# y = distance from the best log-likelihood = (best exact log-likelihood ever
# achieved for that dataset, i.e. the max over every iterate of BOTH rank-1
# curves) - (current exact log-likelihood), on a LOG1P axis. One reference per
# dataset, so both curves are measured against the same target.
#
# y is log1p rather than log10 because the reference is attained BY one of the
# plotted curves, at one of its iterates, so that curve's excess hits exactly 0
# there. Under log10 that point had to be dropped and its neighbours plotted at
# ~1e-3, producing a dramatic-looking plunge whose depth was pure floating-point
# distance-from-itself -- most visible in BBC rank-1 / approximate, and at the
# right edge of pancreas rank-1 / exact. log1p makes 0 representable and
# compresses everything below ~1 toward it, which is also the honest reading: a
# log-likelihood difference under 1 unit is nothing. This is the same argument as
# for the time axis, and it avoids having to invent a best + epsilon offset.
#
# Every panel gets its own x AND y scale; the datasets span very different
# ranges, so compare axis labels across panels, not curve heights.
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
       title = "MOCA K = 25",
       tdiv = 3600, tlab = "Wall-Clock Time (Hours)",   xlog1p = FALSE),
  list(key = "pancreas", tag = "exact_vs_approx_pancreas_K13_c1",
       title = "Pancreas K = 13",
       tdiv = 3600, tlab = "Wall-Clock Time (Hours)",   xlog1p = FALSE),
  list(key = "bbc",      tag = "exact_vs_approx_bbc_K10_c1",
       title = "BBC K = 10",
       tdiv = 1,    tlab = "Wall-Clock Time (Seconds)", xlog1p = TRUE)
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

## same idea on y, but those ranges span up to 13 decades, so thin the grid to
## about five labels
log1p_y_breaks <- function(lims) {
  hi <- max(lims, na.rm = TRUE)
  if (!is.finite(hi) || hi <= 0) return(0)
  top  <- ceiling(log10(max(hi, 1)))
  step <- max(1, ceiling((top + 1) / 5))
  ## start the decades at 10^step, not 10^0: log1p(0) = 0 and log1p(1) = 0.69 sit
  ## on top of each other on an axis that runs to log1p(1e13) = 30, so a break at
  ## 1 just collides with the 0 label
  b    <- c(0, 10^(seq(step, top, by = step)))
  b[b <= hi]
}

sci_lab <- function(x) ifelse(is.na(x), "",
                       ifelse(x == 0, "0", formatC(x, format = "e", digits = 0)))

SCHEME_LABS <- c(random = "Random init", rank1 = "Rank-1 init")
## approx vs exact -- validated colorblind-safe pair (matches the other timing figs)
MCOL <- c("Approximate" = "red", "Exact" = "blue")   # match computational_complexity_fig.R

load_dataset <- function(tag) {
  d <- do.call(rbind, lapply(c("approx", "exact"), function(m) {
    f <- file.path(exp_dir, sprintf("%s_rank1_%s_traces.rds", tag, m))
    if (!file.exists(f)) { message("  missing: ", basename(f)); return(NULL) }
    readRDS(f)
  }))
  d[is.finite(d$exact_loglik), ]
}

## one dataset -> one panel (rank-1 only), no legend
panel_for <- function(ds) {
  d <- load_dataset(ds$tag)
  ## best ever achieved for this dataset: the max over every iterate of both
  ## rank-1 curves. The reference is attained BY one of the plotted curves, at
  ## one of its iterates, which is why one curve per panel touches zero; log1p
  ## renders that as "reached the best" rather than a spurious plunge.
  best     <- max(d$exact_loglik)
  d$excess <- best - d$exact_loglik                # >= 0 by construction
  ## log1p keeps t = 0, so the initialization iterate stays in the plot
  d$hours  <- d$time / ds$tdiv
  d$method <- factor(ifelse(d$method == "approx", "Approximate", "Exact"),
                     levels = c("Approximate", "Exact"))
  g <- ggplot(d, aes(hours, excess, color = method)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = MCOL, name = "Optimization") +
    scale_y_continuous(trans = "log1p", breaks = log1p_y_breaks, labels = sci_lab) +
    labs(x = ds$tlab, y = "Distance from Best Log-Likelihood",
         title = ds$title) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "none")
  if (ds$xlog1p)
    g <- g + scale_x_continuous(trans = "log1p", breaks = log1p_breaks,
                                labels = dec_lab)
  g
}

panels <- lapply(DATASETS, panel_for)

## shared horizontal legend below the panels, left-aligned (as in the
## computational-complexity figure)
legend <- get_legend(panels[[1]] +
                       theme(legend.position = "bottom",
                             legend.justification = "left"))

grid  <- plot_grid(plotlist = panels, nrow = 1, align = "h", axis = "tb",
                   labels = c("A", "B", "C"), label_size = 13)
final <- plot_grid(grid, legend, ncol = 1, rel_heights = c(1, 0.09))

ggsave(out_png, final, width = 12.5, height = 4.3, dpi = 150, bg = "white")
message("Wrote ", normalizePath(out_png, mustWork = FALSE))
