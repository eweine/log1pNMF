# Visualization for note_c0_geometry.tex (the c -> 0+ geometry note in
# the overleaf repo): a worked p = 3 example where the simplex of
# achievable lambda-directions is a drawable triangle.
#
#   A  factor space: the p = 3 structural rows (f_j1, f_j2), with the
#      upper right convex hull boundary; here all three rows are
#      vertices, so the theorem's sparse set consists of the two
#      segments e1--e2 and e2--e3.
#   B  the simplex: the intercept-free theorem allows ONLY the red
#      boundary segments; the intercept-augmented model achieves the
#      whole shaded region cl{softmax(F_s l_s)} -- dense interior, with
#      the red set as its boundary at infinity. Thin curves are the
#      ray paths l_s = t * mu for fixed directions mu, flowing from the
#      uniform vector (t = 0) to their concentration limits.
#
# Example: f_.1 = (3, 2, 0.5), f_.2 = (0.5, 2, 3) across features
# j = 1, 2, 3 (rows in general position; all on the upper right hull).

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))

F1 <- c(3, 2, 0.5)
F2 <- c(0.5, 2, 3)

softmax <- function(e) { e <- e - max(e); exp(e) / sum(exp(e)) }
lam_of  <- function(l1, l2) softmax(l1 * F1 + l2 * F2)

## barycentric projection: e1 -> (0,0), e2 -> (1,0), e3 -> (1/2, sqrt(3)/2)
V <- rbind(c(0, 0), c(1, 0), c(1/2, sqrt(3)/2))
proj <- function(lam) as.numeric(t(V) %*% lam)

## ---- Panel A: factor space with upper right hull ---------------------------
pa_df <- data.frame(x = F1, y = F2, lab = paste0("j = ", 1:3))
hull  <- data.frame(x = c(0.5, 2, 3), y = c(3, 2, 0.5))  # ordered by x
pA <- ggplot(pa_df, aes(x, y)) +
  geom_path(data = hull, colour = "#0072B2", linewidth = 0.9) +
  geom_point(size = 2.6) +
  geom_text(aes(label = lab), nudge_x = 0.16, nudge_y = 0.16, size = 4.2) +
  coord_equal(xlim = c(0, 3.5), ylim = c(0, 3.5)) +
  labs(x = expression(f[j1]), y = expression(f[j2]),
       title = "Structural rows and upper right hull") +
  theme(plot.title = element_text(size = 12.5, hjust = 0.5))

## ---- Panel B: the simplex ---------------------------------------------------
curve_pts <- function(f, tmax = 40, n = 400) {
  t(sapply(seq(0, 1, length.out = n)^2 * tmax, function(t)
    proj(softmax(t * f))))
}
c1 <- curve_pts(F1)                    # center -> e1  (l2 = 0 axis)
c2 <- curve_pts(F2)                    # center -> e3  (l1 = 0 axis)

## achievable region: axis curve to e1, edge e1->e2, edge e2->e3,
## reversed axis curve e3 -> center
poly <- rbind(c1, V[2, ], V[3, ], c2[rev(seq_len(nrow(c2))), ])
poly_df <- data.frame(x = poly[, 1], y = poly[, 2])

tri <- data.frame(x = V[c(1, 2, 3, 1), 1], y = V[c(1, 2, 3, 1), 2])
frontier <- data.frame(x = V[c(1, 2, 3), 1], y = V[c(1, 2, 3), 2])

## ray paths l_s = t * mu: generic directions converge to vertices; the
## two tie-perpendicular directions land inside the segments
tie12 <- c(1.5, 1); tie12 <- tie12 / sqrt(sum(tie12^2))   # perp to row2-row1
tie23 <- c(1, 1.5); tie23 <- tie23 / sqrt(sum(tie23^2))   # perp to row3-row2
rays  <- list(c(1, 0.15), tie12, c(1, 1), tie23, c(0.15, 1))
ray_df <- do.call(rbind, lapply(seq_along(rays), function(i) {
  m <- rays[[i]] / sqrt(sum(rays[[i]]^2))
  pts <- t(sapply(seq(0, 1, length.out = 200)^2 * 40, function(t)
    proj(softmax(t * (m[1] * F1 + m[2] * F2)))))
  data.frame(x = pts[, 1], y = pts[, 2], g = i)
}))

ctr <- proj(rep(1/3, 3))
lab_df <- data.frame(
  x = c(V[1,1] - 0.045, V[2,1] + 0.055, V[3,1], ctr[1] - 0.21),
  y = c(V[1,2] - 0.035, V[2,2] - 0.035, V[3,2] + 0.05, ctr[2] + 0.02),
  lab = c("e[1]", "e[2]", "e[3]", "uniform~(l[s] == 0)"))

pB <- ggplot() +
  geom_polygon(data = poly_df, aes(x, y), fill = "#9BB7BD", alpha = 0.45) +
  geom_path(data = tri, aes(x, y), colour = "grey55", linewidth = 0.5) +
  geom_path(data = ray_df, aes(x, y, group = g), colour = "grey35",
            linewidth = 0.35) +
  geom_path(data = frontier, aes(x, y), colour = "#C2410C",
            linewidth = 1.6) +
  geom_path(data = data.frame(x = c1[, 1], y = c1[, 2]), aes(x, y),
            colour = "#0E5C6B", linewidth = 0.8) +
  geom_path(data = data.frame(x = c2[, 1], y = c2[, 2]), aes(x, y),
            colour = "#0E5C6B", linewidth = 0.8) +
  geom_point(aes(x = ctr[1], y = ctr[2]), size = 2.2) +
  geom_text(data = lab_df, aes(x, y, label = lab), parse = TRUE,
            size = 4.2) +
  coord_equal(xlim = c(-0.1, 1.1), ylim = c(-0.08, 0.95)) +
  labs(title = "Achievable directions of λ on the simplex") +
  theme_void() +
  theme(plot.title = element_text(size = 12.5, hjust = 0.5))

g <- plot_grid(pA, pB, nrow = 1, labels = c("A", "B"), label_size = 13,
               rel_widths = c(0.85, 1))
ggsave("../images/note_c0_geometry_fig.png", g, width = 9.6, height = 4.4,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_geometry_fig.png")
