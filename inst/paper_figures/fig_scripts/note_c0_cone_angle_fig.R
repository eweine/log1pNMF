# Diagram for Proposition 3 ("Area supremum, p = 3") of
# note_c0_geometry.tex: the opening angle of the loading cone drives
# the area of the achievable set.
#
# Two configurations, one per row:
#   A  small area: rows near an INCREASING line, so the centered
#      columns g1, g2 are nearly parallel and the loading cone is a
#      needle (theta ~ 10 deg); F = cbind(c(1,2,2.1), c(1,2.2,2.1)).
#   B  large area: the near-maximizing family of Proposition 3 at
#      s = 0.05 -- rows near a DECREASING line, centered columns
#      nearly antipodal, cone nearly a half-plane (theta ~ 170 deg);
#      rows (2.05, 0.05), (0.05, 2.05), (1.20, 1.20).
#
# Columns: factor space (rows of F) | loading cone in the centered
# plane (the wedge cone(g1, g2), with its opening angle) | achievable
# set on the simplex (with exact area from the rejection sampler).

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 11))

V <- rbind(c(0, 0), c(1, 0), c(1/2, sqrt(3)/2))
proj_mat <- function(P) P %*% V
softmax_rows <- function(E) {
  E <- E - apply(E, 1, max)
  W <- exp(E); W / rowSums(W)
}
region_cloud <- function(F, Tmax, n = 420) {
  g  <- seq(0, 1, length.out = n)^1.6 * Tmax
  ll <- as.matrix(expand.grid(l1 = g, l2 = g))
  proj_mat(softmax_rows(ll %*% t(F)))
}
area_frac <- function(F, n = 6e5) {
  set.seed(1)
  E <- matrix(rexp(3 * n), ncol = 3); L <- E / rowSums(E)
  eta <- log(L); eta <- eta - rowMeans(eta)
  B <- cbind(F[, 1] - mean(F[, 1]), F[, 2] - mean(F[, 2]))
  cf <- t(qr.solve(B, t(eta)))
  mean(cf[, 1] >= -1e-9 & cf[, 2] >= -1e-9)
}

## orthonormal basis of the centered plane, for drawing g1, g2 in 2-D
b1 <- c(1, -1, 0) / sqrt(2)
b2 <- c(1, 1, -2) / sqrt(6)
to2d <- function(g) c(sum(g * b1), sum(g * b2))

panel_factor <- function(F, line_pts, title) {
  d <- data.frame(x = F[, 1], y = F[, 2], lab = paste0("j=", 1:3))
  lim <- c(0, max(F) * 1.25)
  ggplot(d, aes(x, y)) +
    geom_segment(aes(x = line_pts[1], y = line_pts[2],
                     xend = line_pts[3], yend = line_pts[4]),
                 colour = "grey55", linetype = "dashed") +
    geom_point(size = 2.4) +
    geom_text(aes(label = lab), nudge_x = diff(lim) * 0.07,
              nudge_y = diff(lim) * 0.05, size = 3.4) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = expression(f[j1]), y = expression(f[j2]), title = title) +
    theme(plot.title = element_text(size = 11, hjust = 0.5))
}

panel_cone <- function(F, title) {
  g1 <- to2d(F[, 1] - mean(F[, 1]))
  g2 <- to2d(F[, 2] - mean(F[, 2]))
  a1 <- atan2(g1[2], g1[1]); a2 <- atan2(g2[2], g2[1])
  th <- acos(sum(g1 * g2) / sqrt(sum(g1^2) * sum(g2^2)))
  ## wedge (short way from a1 to a2, angle th < pi)
  dlt <- a2 - a1
  if (dlt >  pi) dlt <- dlt - 2 * pi
  if (dlt < -pi) dlt <- dlt + 2 * pi
  aa <- a1 + seq(0, 1, length.out = 120) * dlt
  wedge <- data.frame(x = c(0, cos(aa)), y = c(0, sin(aa)))
  arc <- data.frame(x = 0.32 * cos(aa), y = 0.32 * sin(aa))
  u1 <- g1 / sqrt(sum(g1^2)); u2 <- g2 / sqrt(sum(g2^2))
  mid <- a1 + dlt / 2
  ggplot() +
    geom_polygon(data = wedge, aes(x, y), fill = "#9BB7BD", alpha = 0.45) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_path(data = arc, aes(x, y), colour = "#C2410C", linewidth = 0.7) +
    geom_segment(aes(x = 0, y = 0, xend = u1[1], yend = u1[2]),
                 colour = "#0E5C6B", linewidth = 0.9,
                 arrow = arrow(length = unit(0.09, "in"))) +
    geom_segment(aes(x = 0, y = 0, xend = u2[1], yend = u2[2]),
                 colour = "#0E5C6B", linewidth = 0.9,
                 arrow = arrow(length = unit(0.09, "in"))) +
    annotate("text", x = 1.13 * u1[1], y = 1.13 * u1[2],
             label = "g[1]", parse = TRUE, size = 3.6, colour = "#0E5C6B") +
    annotate("text", x = 1.13 * u2[1], y = 1.13 * u2[2],
             label = "g[2]", parse = TRUE, size = 3.6, colour = "#0E5C6B") +
    { # place the angle label inside a wide wedge, beside a narrow one
      pos <- if (th > pi/3) 0.5 * c(cos(mid), sin(mid)) else {
        perp <- c(-sin(mid), cos(mid))
        if (perp[2] < 0) perp <- -perp
        0.42 * c(cos(mid), sin(mid)) + 0.34 * perp
      }
      annotate("text", x = pos[1], y = pos[2],
               label = sprintf("theta %%~~%% '%.0f'*degree",
                               th * 180 / pi),
               parse = TRUE, size = 3.6, colour = "#C2410C")
    } +
    coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
    labs(title = title) + theme_void() +
    theme(plot.title = element_text(size = 11, hjust = 0.5))
}

panel_simplex <- function(F, title) {
  pts <- region_cloud(F, Tmax = 300)
  tri <- data.frame(x = V[c(1, 2, 3, 1), 1], y = V[c(1, 2, 3, 1), 2])
  u <- proj_mat(matrix(1/3, 1, 3))
  ggplot() +
    geom_path(data = tri, aes(x, y), colour = "grey55", linewidth = 0.4) +
    geom_point(data = data.frame(x = pts[, 1], y = pts[, 2]),
               aes(x, y), colour = "#9BB7BD", size = 0.2, alpha = 0.3) +
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
}

FA <- cbind(c(1, 2, 2.1), c(1, 2.2, 2.1))         # needle
FB <- cbind(c(2.05, 0.05, 1.20), c(0.05, 2.05, 1.20))  # near-maximizer, s = 0.05

rowA <- plot_grid(
  panel_factor(FA, c(0.9, 0.9, 2.2, 2.35), "rows near an increasing line"),
  panel_cone(FA, "loading cone: needle"),
  panel_simplex(FA, sprintf("area = %.3f |Δ|", area_frac(FA))),
  nrow = 1, rel_widths = c(0.9, 1, 1))
rowB <- plot_grid(
  panel_factor(FB, c(2.15, -0.05, -0.05, 2.15), "rows near a decreasing line"),
  panel_cone(FB, "loading cone: near half-plane"),
  panel_simplex(FB, sprintf("area = %.3f |Δ|", area_frac(FB))),
  nrow = 1, rel_widths = c(0.9, 1, 1))
g <- plot_grid(rowA, rowB, ncol = 1, labels = c("A", "B"), label_size = 12)
ggsave("../images/note_c0_cone_angle_fig.png", g, width = 9.6, height = 6.6,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_cone_angle_fig.png")
