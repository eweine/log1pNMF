# Structure plots of the loadings from one c-grid cell: true L above fitted L.
#
#   Rscript peek_c_grid_L.R <seed> <c_true> <c_fit> [init=] [rownorm=] [out=]
#
#   Rscript peek_c_grid_L.R 1 1 1                      # well-specified
#   Rscript peek_c_grid_L.R 1 1 Inf                    # same data, fit at c = Inf
#   Rscript peek_c_grid_L.R 3 0.001 1000 init=random
#   Rscript peek_c_grid_L.R 1 1 Inf rownorm=1 out=x.png
#
# c_true / c_fit take the grid values 0.001 0.01 0.1 1 10 100 1000 Inf.
# init is "rank1" (default) or "random".
# Point C_GRID_OUT at the cell directory if it is not the cluster default.
#
# NORMALIZATION
#
# The column step is the paper convention (normalize_bars): divide each loading
# column by its maximum so every factor peaks at one. Without it a factor with
# systematically small loadings is invisible, because each cell's bar is
# dominated by whichever factor carries the largest raw values.
#
# rownorm=0 (default) reproduces normalized_structure_plot exactly. Note this
# does NOT give memberships: fastTopics::structure_plot row-normalizes only when
# handed a fit object (it runs poisson2multinom); given a plain matrix it plots
# the values as they are. So bars vary in total height -- verified directly, a
# matrix with row sums 1.4, 0.76, 1.9 plots at exactly those heights. The true L
# here happens to land near 1.0 only because its rows are Dirichlet draws and
# each column max is near one; the fitted L reaches about 1.2.
#
# rownorm=1 divides each row by its sum after the column step, so both panels
# are genuine membership proportions on a common 0-1 scale. Use it when
# comparing the two panels; use the default when matching the paper figures.
#
# Both panels use the SAME cell ordering (by the true dominant factor, then by
# that factor's weight) and the same colours, and the fitted factors are
# permuted to the truth by the rule score_fit uses -- the column permutation
# maximizing total correlation with the true L. Without the permutation the two
# panels differ by an arbitrary relabelling and every fit looks wrong.

suppressMessages({
  library(fastTopics)
  library(log1pNMF)
  library(ggplot2)
  library(cowplot)
})
source("c_grid_sim_common.R")

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3)
  stop("usage: Rscript peek_c_grid_L.R <seed> <c_true> <c_fit> [init=] [rownorm=] [out=]")

kv <- a[grepl("=", a)]
opt <- setNames(sub("^[^=]*=", "", kv), sub("=.*$", "", kv))
getopt <- function(k, default) if (k %in% names(opt)) opt[[k]] else default

seed   <- as.integer(a[1])
c_true <- as.numeric(a[2])
c_fit  <- as.numeric(a[3])
init   <- getopt("init", "rank1")
rownrm <- as.integer(getopt("rownorm", "0")) == 1L
fmtc   <- function(x) format(x, scientific = FALSE, trim = TRUE)
out    <- getopt("out",
  sprintf("peek_L_seed%02d_ctrue%s_cfit%s_%s%s.png", seed, fmtc(c_true), fmtc(c_fit),
          init, if (rownrm) "_rownorm" else ""))

spec <- list(seed = seed, c_true = c_true, c_fit = c_fit)
f    <- file.path(OUT_DIR, paste0(tag_of(spec), ".rds"))
if (!file.exists(f)) stop("no such cell: ", f)
x <- readRDS(f)
if (is.null(x$fits[[init]]))
  stop("no init '", init, "' in that cell; have: ", paste(names(x$fits), collapse = ", "))

Lt <- x$L_true
Lf <- x$fits[[init]]$LL[, match_factors(Lt, x$fits[[init]]$LL), drop = FALSE]

row <- x$rows[x$rows$init == init, ]
cat(sprintf("seed %d | c_true %s | c_fit %s | init %s\n",
            seed, fmtc(c_true), fmtc(c_fit), init))
cat(sprintf("  L_cor_mean %.4f  F_cor_mean %.4f  loglik %.1f  n_iter %d  converged %s\n",
            row$L_cor_mean, row$F_cor_mean, row$loglik, row$n_iter, row$converged))
cat("  per-factor correlation with the truth: ",
    paste(sprintf("k%d %.3f", seq_len(ncol(Lt)),
                  diag(suppressWarnings(cor(Lt, Lf)))), collapse = "   "), "\n", sep = "")

## shared ordering: by the true dominant factor, then that factor's weight
dom <- max.col(Lt)
ord <- order(dom, -Lt[cbind(seq_along(dom), dom)])
grp <- factor(dom[ord])
KK  <- ncol(Lt)
TOP <- paste0("k", seq_len(KK))
## validated colourblind-safe triple (all six checks pass on a light surface)
FCOL <- setNames(c("#0E8FA8", "#C2410C", "#6D3FC4")[seq_len(KK)], TOP)

sp <- function(M, title, show_x = FALSE) {
  M <- pmax(M[ord, , drop = FALSE], 0)
  cm <- apply(M, 2, max)
  if (any(cm <= 0)) stop("a loading column is entirely zero; cannot normalize to max 1")
  Mn <- sweep(M, 2, cm, "/")                       # normalize_bars: column max = 1
  if (rownrm) Mn <- Mn / pmax(rowSums(Mn), .Machine$double.eps)
  colnames(Mn) <- TOP
  g <- structure_plot(Mn, loadings_order = seq_len(nrow(Mn)), grouping = grp,
                      gap = 6, topics = TOP, colors = FCOL) +
    labs(y = if (rownrm) "membership" else "scaled loading", title = title) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.position = "none")
  if (show_x) g + labs(x = "cell, grouped by true dominant factor")
  else        g + theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank()) + labs(x = NULL)
}

leg <- get_legend(
  ggplot(data.frame(x = 1, f = factor(TOP, levels = TOP)), aes(x, 1, fill = f)) +
    geom_col() + scale_fill_manual(values = FCOL, name = "factor") +
    theme_cowplot() + theme(legend.position = "right"))

body <- plot_grid(sp(Lt, "true L"),
                  sp(Lf, sprintf("fitted L   (c_fit = %s, %s init)", fmtc(c_fit), init),
                     show_x = TRUE),
                  ncol = 1, align = "v", axis = "lr")

ttl <- ggdraw() +
  draw_label(sprintf("seed %d, simulated at c = %s   |   mean loading correlation %.3f",
                     seed, fmtc(c_true), row$L_cor_mean),
             x = 0.01, hjust = 0, size = 12)

final <- plot_grid(ttl,
                   plot_grid(body, leg, ncol = 2, rel_widths = c(1, 0.11)),
                   ncol = 1, rel_heights = c(0.07, 1))

ggsave(out, final, width = 10, height = 6, dpi = 150, bg = "white")
message("Wrote ", out)
