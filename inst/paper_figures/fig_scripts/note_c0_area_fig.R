# Figures for the "M vs area" section of note_c0_geometry.tex.
#
# Figure 1 (note_c0_area_fig.png): three configurations of the structural
# factors (rows: factor space with upper-right hull; simplex: T = "inf"
# point cloud of the achievable set), chosen to show that the number of
# hull vertices M and the AREA of the achievable set agree at the
# extremes but decouple in between:
#   A  M = 1, fat lens:  f1 = (4,2,1),  f2 = (4,1,2)  (row 1 dominates)
#   B  M = 2, thin:      f1 = (1,2,2.1), f2 = (1,2.2,2.1)
#      (rows nearly on an increasing line, one trade-off pair)
#   C  M = 3, wide:      f1 = (3,2,1),  f2 = (1,2.05,3)
#      (rows on a decreasing line + epsilon; centered columns nearly
#      opposite, so the loading cone is wide and the region fills)
#
# Figure 2 (note_c0_area_fill_fig.png): config C at increasing loading
# budgets T, showing that the full region is approached only at loadings
# of order 1/epsilon: at bounded T the region is a thin band even though
# M = 3 and the T = inf region is wide.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 11))

V <- rbind(c(0, 0), c(1, 0), c(1/2, sqrt(3)/2))
proj_mat <- function(P) P %*% V            # rows of P = simplex points
softmax_rows <- function(E) {
  E <- E - apply(E, 1, max)
  W <- exp(E); W / rowSums(W)
}

region_cloud <- function(F, Tmax, n = 420) {
  g  <- seq(0, 1, length.out = n)^1.6 * Tmax
  ll <- as.matrix(expand.grid(l1 = g, l2 = g))
  proj_mat(softmax_rows(ll %*% t(F)))
}

frontier_verts <- function(F) {
  th  <- seq(0, pi/2, length.out = 4001)
  amx <- sapply(th, function(a) which.max(cos(a) * F[, 1] + sin(a) * F[, 2]))
  rle(amx)$values
}

panel_factor <- function(F, title) {
  d  <- data.frame(x = F[, 1], y = F[, 2], lab = paste0("j=", 1:3))
  vs <- frontier_verts(F)
  lim <- c(0, max(F) * 1.25)
  g <- ggplot(d, aes(x, y)) + geom_point(size = 2.2)
  if (length(vs) >= 2)
    g <- g + geom_path(data = d[vs, ], colour = "#C2410C", linewidth = 0.8)
  g + geom_point(data = d[unique(vs), ], colour = "#C2410C", size = 2.6,
                 shape = 15) +
    geom_text(aes(label = lab), nudge_x = diff(lim) * 0.06,
              nudge_y = diff(lim) * 0.05, size = 3.4) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = expression(f[j1]), y = expression(f[j2]), title = title) +
    theme(plot.title = element_text(size = 11, hjust = 0.5))
}

panel_simplex <- function(F, Tmax, title, cloud_col = "#9BB7BD") {
  pts <- region_cloud(F, Tmax)
  tri <- data.frame(x = V[c(1, 2, 3, 1), 1], y = V[c(1, 2, 3, 1), 2])
  vs  <- frontier_verts(F)
  u   <- proj_mat(matrix(1/3, 1, 3))
  g <- ggplot() +
    geom_path(data = tri, aes(x, y), colour = "grey55", linewidth = 0.4) +
    geom_point(data = data.frame(x = pts[, 1], y = pts[, 2]),
               aes(x, y), colour = cloud_col, size = 0.2, alpha = 0.3)
  if (length(vs) >= 2) {
    fr <- data.frame(x = V[vs, 1], y = V[vs, 2])
    g <- g + geom_path(data = fr, aes(x, y), colour = "#C2410C",
                       linewidth = 1.3)
  }
  g <- g + geom_point(data = data.frame(x = V[unique(vs), 1],
                                        y = V[unique(vs), 2]),
                      aes(x, y), colour = "#C2410C", shape = 15, size = 2) +
    geom_point(aes(x = u[1], y = u[2]), size = 1.6) +
    annotate("text", x = V[1,1] - 0.05, y = V[1,2] - 0.03,
             label = "e[1]", parse = TRUE, size = 3.2) +
    annotate("text", x = V[2,1] + 0.05, y = V[2,2] - 0.03,
             label = "e[2]", parse = TRUE, size = 3.2) +
    annotate("text", x = V[3,1], y = V[3,2] + 0.04,
             label = "e[3]", parse = TRUE, size = 3.2) +
    coord_equal(xlim = c(-0.09, 1.09), ylim = c(-0.07, 0.94)) +
    labs(title = title) + theme_void() +
    theme(plot.title = element_text(size = 11, hjust = 0.5))
  g
}

CONFIGS <- list(
  list(F = cbind(c(4, 2, 1), c(4, 1, 2)),
       lab = "M = 1, wide cone (fat lens)"),
  list(F = cbind(c(1, 2, 2.1), c(1, 2.2, 2.1)),
       lab = "M = 2, narrow cone (thin)"),
  list(F = cbind(c(3, 2, 1), c(1, 2.05, 3)),
       lab = "M = 3, wide cone (fills)")
)

rows <- lapply(CONFIGS, function(cf) {
  plot_grid(panel_factor(cf$F, cf$lab),
            panel_simplex(cf$F, Tmax = 300, title = "achievable set (T = ∞)"),
            nrow = 1, rel_widths = c(0.8, 1))
})
g1 <- plot_grid(plotlist = rows, ncol = 1, labels = c("A", "B", "C"),
                label_size = 12)
ggsave("../images/note_c0_area_fig.png", g1, width = 7.2, height = 9.6,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_area_fig.png")

## Figure 2: config C at increasing budgets
FC <- CONFIGS[[3]]$F
panels <- lapply(c(5, 20, 300), function(Tm)
  panel_simplex(FC, Tmax = Tm,
                title = ifelse(Tm == 300, "T = ∞ (T = 300)",
                               paste0("T = ", Tm))))
g2 <- plot_grid(plotlist = panels, nrow = 1)
ggsave("../images/note_c0_area_fill_fig.png", g2, width = 10.2, height = 3.4,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_area_fill_fig.png")

## Figure 3: the two dominance counterexamples, pairwise
## (note_c0_counterex_fig.png). Top pair: an M = 1 configuration whose
## area exceeds an M = 2 configuration's. Bottom pair: an M = 2
## configuration whose area exceeds an M = 3 configuration's (found by
## random search; areas computed exactly via the softmax-diffeomorphism
## rejection sampler).
area_frac <- function(F, n = 6e5) {
  set.seed(1)
  E <- matrix(rexp(3 * n), ncol = 3); L <- E / rowSums(E)
  eta <- log(L); eta <- eta - rowMeans(eta)
  B <- cbind(F[, 1] - mean(F[, 1]), F[, 2] - mean(F[, 2]))
  cf <- t(qr.solve(B, t(eta)))
  mean(cf[, 1] >= -1e-9 & cf[, 2] >= -1e-9)
}

CE <- list(
  list(F = cbind(c(4, 2, 1), c(4, 1, 2))),                    # M = 1
  list(F = cbind(c(1, 2, 2.1), c(1, 2.2, 2.1))),              # M = 2 thin
  list(F = cbind(c(0.12, 0.19, 2.50), c(3.47, 3.48, 0.09))),  # M = 2 wide
  list(F = cbind(c(3.07, 2.08, 3.05), c(1.29, 3.96, 3.89)))   # M = 3
)
ce_panels <- lapply(CE, function(cf) {
  M <- length(unique(frontier_verts(cf$F)))
  a <- area_frac(cf$F)
  plot_grid(panel_factor(cf$F, "factor space"),
            panel_simplex(cf$F, Tmax = 300,
                          title = sprintf("M = %d:  area = %.3f |Δ|", M, a)),
            nrow = 1, rel_widths = c(0.72, 1))
})
g3 <- plot_grid(plotlist = ce_panels, nrow = 2, labels = c("A", "", "B", ""),
                label_size = 12)
ggsave("../images/note_c0_counterex_fig.png", g3, width = 12.4, height = 6.8,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_counterex_fig.png")
