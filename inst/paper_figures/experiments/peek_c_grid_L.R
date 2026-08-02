# Look at the loadings from one c-grid cell: true L beside the fitted L, with
# the fitted factors permuted to match the truth.
#
#   Rscript peek_c_grid_L.R <seed> <c_true> <c_fit> [init] [out.png]
#
#   Rscript peek_c_grid_L.R 1 1 1                 # well-specified, rank-1 init
#   Rscript peek_c_grid_L.R 1 1 Inf random        # same data, fit at c = Inf
#   Rscript peek_c_grid_L.R 3 0.001 1000
#
# c_true / c_fit take the grid values 0.001 0.01 0.1 1 10 100 1000 Inf.
# init is "rank1" (default) or "random".
#
# Cells are sorted by their TRUE dominant factor, then by that factor's weight,
# and the same ordering is used for both panels -- so a clean fit shows three
# blocks matching the truth, and a poor one shows the blocks smearing together.
#
# The fitted factors are permuted by the same rule score_fit uses: the column
# permutation maximizing total correlation with the true L. Without it the
# panels would differ by an arbitrary relabelling and look worse than they are.

suppressMessages({ library(ggplot2); library(cowplot) })
theme_set(theme_cowplot(font_size = 11))
source("c_grid_sim_common.R")

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3)
  stop("usage: Rscript peek_c_grid_L.R <seed> <c_true> <c_fit> [init] [out.png]")

seed   <- as.integer(a[1])
c_true <- as.numeric(a[2])
c_fit  <- as.numeric(a[3])
init   <- if (length(a) >= 4) a[4] else "rank1"
out    <- if (length(a) >= 5) a[5] else
  sprintf("peek_L_seed%02d_ctrue%s_cfit%s_%s.png", seed,
          format(c_true, scientific = FALSE, trim = TRUE),
          format(c_fit,  scientific = FALSE, trim = TRUE), init)

spec <- list(seed = seed, c_true = c_true, c_fit = c_fit)
f    <- file.path(OUT_DIR, paste0(tag_of(spec), ".rds"))
if (!file.exists(f)) stop("no such cell: ", f)
x <- readRDS(f)
if (is.null(x$fits[[init]])) stop("no init '", init, "' in that cell")

Lt <- x$L_true
Lf <- x$fits[[init]]$LL
pp <- match_factors(Lt, Lf)          # align fitted columns to the truth
Lf <- Lf[, pp, drop = FALSE]

row <- x$rows[x$rows$init == init, ]
cat(sprintf("seed %d | c_true %s | c_fit %s | init %s\n", seed,
            format(c_true), format(c_fit), init))
cat(sprintf("  L_cor_mean %.4f  F_cor_mean %.4f  loglik %.1f  n_iter %d  converged %s\n",
            row$L_cor_mean, row$F_cor_mean, row$loglik, row$n_iter, row$converged))
cat(sprintf("  factor permutation applied to the fit: %s\n", paste(pp, collapse = "-")))
cat("\n  per-factor correlation with the truth: ",
    paste(sprintf("k%d %.3f", seq_len(ncol(Lt)),
                  diag(suppressWarnings(cor(Lt, Lf)))), collapse = "   "), "\n", sep = "")

## rows ordered by the TRUE dominant factor, then that factor's weight
dom <- max.col(Lt)
ord <- order(dom, -Lt[cbind(seq_len(nrow(Lt)), dom)])

long <- function(M, lab) {
  M <- M[ord, , drop = FALSE]
  M <- sweep(M, 1, pmax(rowSums(M), .Machine$double.eps), "/")   # row proportions
  data.frame(cell = rep(seq_len(nrow(M)), ncol(M)),
             factor = factor(rep(seq_len(ncol(M)), each = nrow(M))),
             w = as.vector(M), panel = lab)
}
d <- rbind(long(Lt, "true L"),
           long(Lf, sprintf("fitted L  (c_fit = %s, %s)", format(c_fit), init)))
d$panel <- factor(d$panel, levels = unique(d$panel))

## three factors -> a validated colourblind-safe triple
FCOL <- c("1" = "#0E8FA8", "2" = "#C2410C", "3" = "#6D3FC4")

p <- ggplot(d, aes(cell, w, fill = factor)) +
  geom_col(width = 1) +
  facet_wrap(~ panel, ncol = 1) +
  scale_fill_manual(values = FCOL, name = "factor") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "cell (ordered by true dominant factor)", y = "membership",
       title = sprintf("Loadings, seed %d, simulated at c = %s",
                       seed, format(c_true)),
       subtitle = sprintf("mean correlation with the truth: %.3f", row$L_cor_mean)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 12, face = "plain"),
        plot.subtitle = element_text(size = 10, colour = "grey35"))

ggsave(out, p, width = 9, height = 5.5, dpi = 150, bg = "white")
message("Wrote ", out)
