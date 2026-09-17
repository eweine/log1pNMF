# Explanatory diagrams for the proofs in note_c0_geometry.tex.
#
#   Figure 1 (note_c0_exhaust_fig.png): the half-plane exhaustion
#     (Lemma 2). In the plane V with coordinates (a, b) along the
#     orthonormal basis (w, u): the cone of g1 = w + eps*u and
#     g2 = -w + eps*u is the wedge {b >= eps*|a|}, which opens toward
#     the half-plane H = {b >= 0} as eps decreases.
#
#   Figure 2 (note_c0_join_fig.png): the join surface of Lemma 4 at
#     p = 4 (the simplex is a solid tetrahedron, drawn in projection).
#     The log-line gamma(s) lies in the bottom facet {lambda_4 = 0};
#     the surface is the union of the rulings (segments) from the
#     apex e4 to the curve; the waist curve b = 0 (t = t0(s)) passes
#     through the uniform vector and splits the surface into the apex
#     half (b >= 0, shaded) and the facet half.
#
#   Figure 3 (note_c0_crawl_fig.png): the crawl (Section 5), p = 4.
#     Left: the log-line in the facet triangle for d_delta =
#     (2, -1+delta, -1-delta) at decreasing delta: the curve runs from
#     e1 through the centroid toward mid(e2,e3) and then crawls along
#     the edge toward e3; in the limit it is the symmetric arc plus
#     the half-edge. Right: the limit surface in the tetrahedron: the
#     join over the arc plus the gained flat triangle
#     conv(e4, mid(e2,e3), e3) -- half a 2-face, area sqrt(3)/4.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 11))

TEAL  <- "#0E5C6B"; CLOUD <- "#9BB7BD"; ORANGE <- "#C2410C"
GREY  <- "grey55"

softmax_rows <- function(E) {
  E <- E - apply(E, 1, max)
  W <- exp(E); W / rowSums(W)
}

## ---------- Figure 1: exhaustion wedges -------------------------------------
wedge_poly <- function(eps, A = 1.25, B = 1.25) {
  data.frame(a = c(-A, 0, A, A, -A), b = c(eps * A, 0, eps * A, B, B),
             eps = eps)
}
eps_vals <- c(0.95, 0.45, 0.15)
polys <- do.call(rbind, lapply(eps_vals, wedge_poly))
polys$eps <- factor(polys$eps, levels = sort(eps_vals, decreasing = TRUE))
gvec <- function(eps) { g <- c(1, eps); 1.05 * g / sqrt(sum(g^2)) }
g1 <- gvec(0.45); g2 <- c(-1, 1) * gvec(0.45)
p1 <- ggplot() +
  geom_polygon(data = polys, aes(a, b, group = eps, alpha = eps),
               fill = CLOUD) +
  scale_alpha_manual(values = c(0.25, 0.25, 0.25), guide = "none") +
  geom_hline(yintercept = 0, colour = TEAL, linewidth = 0.9) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_segment(aes(x = 0, y = 0, xend = g1[1], yend = g1[2]),
               colour = TEAL, linewidth = 0.9,
               arrow = arrow(length = unit(0.09, "in"))) +
  geom_segment(aes(x = 0, y = 0, xend = g2[1], yend = g2[2]),
               colour = TEAL, linewidth = 0.9,
               arrow = arrow(length = unit(0.09, "in"))) +
  annotate("text", x = g1[1] + 0.14, y = g1[2] + 0.02,
           label = "g[1] == w + epsilon*minute*u", parse = TRUE,
           size = 3.4, colour = TEAL) +
  annotate("text", x = g2[1] - 0.14, y = g2[2] + 0.02,
           label = "g[2] == -w + epsilon*minute*u", parse = TRUE,
           size = 3.4, colour = TEAL) +
  annotate("text", x = -0.62, y = 1.2,
           label = "cone(g[1], g[2]) == group('{', b >= epsilon*minute*'|'*a*'|', '}')",
           parse = TRUE, size = 3.6) +
  annotate("text", x = 0.3, y = 0.62,
           label = "epsilon*minute %down% 0",
           parse = TRUE, size = 3.8, colour = ORANGE) +
  annotate("text", x = -1.02, y = -0.09,
           label = "boundary~of~H*':'~the~line~b == 0", parse = TRUE,
           size = 3.3, colour = TEAL) +
  annotate("text", x = 1.08, y = -0.08, label = "a  (along w)", size = 3.3) +
  annotate("text", x = 0.08, y = 0.9, label = "b  (along u)", size = 3.3,
           angle = 90, colour = "grey40") +
  coord_equal(xlim = c(-1.3, 1.3), ylim = c(-0.18, 1.3)) +
  theme_void()
ggsave("../images/note_c0_exhaust_fig.png", p1, width = 6.4, height = 4.0,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_exhaust_fig.png")

## ---------- 3-D tetrahedron machinery (p = 4) -------------------------------
## regular tetrahedron vertices in R^3 (rows = e1..e4), then an
## orthographic projection to 2-D chosen so e4 is at the top and the
## facet {1,2,3} is the floor.
T3 <- rbind(c( 1,  1, -1),
            c(-1, -1, -1),
            c( 1, -1,  1),
            c(-1,  1,  1)) / sqrt(3)
## rotate so e4 points up: align e4 with +z
alignz <- function(v) {
  v <- v / sqrt(sum(v^2)); z <- c(0, 0, 1)
  ax <- c(v[2]*z[3]-v[3]*z[2], v[3]*z[1]-v[1]*z[3], v[1]*z[2]-v[2]*z[1])
  s <- sqrt(sum(ax^2)); cthet <- sum(v*z)
  if (s < 1e-12) return(diag(3))
  ax <- ax / s; K <- rbind(c(0,-ax[3],ax[2]), c(ax[3],0,-ax[1]),
                           c(-ax[2],ax[1],0))
  diag(3) + sqrt(1-cthet^2)*K + (1-cthet)*(K %*% K)
}
R1 <- alignz(T3[4, ] - colMeans(T3))
T3 <- T3 %*% t(R1)
az <- 100 * pi/180
Rz <- rbind(c(cos(az), -sin(az), 0), c(sin(az), cos(az), 0), c(0, 0, 1))
T3 <- T3 %*% t(Rz)
el <- 14 * pi/180   # camera elevation
v1 <- c(1, 0, 0)
v2 <- c(0, sin(el), cos(el))
proj4 <- function(L) {               # rows: lambda (4 cols, sum 1)
  X <- L %*% T3
  data.frame(x = as.numeric(X %*% v1), y = as.numeric(X %*% v2))
}
E4 <- diag(4)
tet_edges <- function() {
  idx <- combn(4, 2)
  do.call(rbind, lapply(seq_len(ncol(idx)), function(k) {
    a <- proj4(E4[idx[1, k], , drop = FALSE])
    b <- proj4(E4[idx[2, k], , drop = FALSE])
    data.frame(x = a$x, y = a$y, xend = b$x, yend = b$y)
  }))
}
vlab <- data.frame(proj4(E4),
                   lab = c("e[1]", "e[2]", "e[3]", "e[4]"))
vlab$y <- vlab$y + c(-0.09, -0.09, -0.09, 0.09)

gamma_curve <- function(d, s) {      # d: length 3 (facet), rows of lambda
  G <- softmax_rows(outer(s, d))
  cbind(G, 0)
}
join_pts <- function(d, s, tt) {     # full (s, t) grid of surface points
  Gam <- gamma_curve(d, s)
  do.call(rbind, lapply(tt, function(t)
    cbind((1 - t) * Gam[, 1:3], t)))
}

## ---------- Figure 2: the join surface --------------------------------------
d0 <- c(2, -1, -1)
sg <- tan(seq(-1, 1, length.out = 241) * (pi/2 - 0.02)) * 3
tt <- seq(0.02, 0.98, length.out = 60)
Gam <- gamma_curve(d0, sg)
Z   <- rowSums(exp(outer(sg, d0)))
t0  <- 1 / (1 + Z)
## surface cloud split by half
cloud <- do.call(rbind, lapply(seq_along(tt), function(i) {
  P <- cbind((1 - tt[i]) * Gam[, 1:3], tt[i])
  data.frame(proj4(P), apex = tt[i] >= t0)
}))
waist <- data.frame(proj4(cbind((1 - t0) * Gam[, 1:3], t0)))
curveF <- data.frame(proj4(Gam))
rul_s <- sg[round(seq(6, length(sg) - 5, length.out = 13))]
ruls <- do.call(rbind, lapply(rul_s, function(s) {
  g <- gamma_curve(d0, s); a <- proj4(g); b <- proj4(E4[4, , drop = FALSE])
  data.frame(x = a$x, y = a$y, xend = b$x, yend = b$y)
}))
u4 <- proj4(matrix(1/4, 1, 4))
p2 <- ggplot() +
  geom_segment(data = tet_edges(), aes(x, y, xend = xend, yend = yend),
               colour = "grey75", linewidth = 0.35) +
  geom_point(data = subset(cloud, !apex), aes(x, y), colour = "grey80",
             size = 0.25, alpha = 0.35) +
  geom_point(data = subset(cloud, apex), aes(x, y), colour = CLOUD,
             size = 0.3, alpha = 0.5) +
  geom_segment(data = ruls, aes(x, y, xend = xend, yend = yend),
               colour = "grey45", linewidth = 0.25) +
  geom_path(data = curveF, aes(x, y), colour = TEAL, linewidth = 1.1) +
  geom_path(data = waist, aes(x, y), colour = ORANGE, linewidth = 0.9) +
  geom_point(aes(x = u4$x, y = u4$y), size = 1.8) +
  geom_text(data = vlab, aes(x, y, label = lab), parse = TRUE, size = 3.8) +
  annotate("text", x = u4$x - 0.33, y = u4$y + 0.02,
           label = "uniform~(s == 0*','~b == 0)", parse = TRUE, size = 3.2) +
  annotate("text", x = -0.1, y = -0.8,
           label = "log*'-'*line~gamma(s)~'in the facet'", parse = TRUE,
           size = 3.4, colour = TEAL) +
  annotate("text", x = -0.78, y = -0.02,
           label = "waist~b == 0*':'~t == t[0](s)", parse = TRUE,
           size = 3.4, colour = ORANGE) +
  annotate("text", x = 0.78, y = 0.42,
           label = "apex~half~(b >= 0)", parse = TRUE, size = 3.4,
           colour = "#5F8790") +
  coord_equal(xlim = c(-1.2, 1.3)) + theme_void()
ggsave("../images/note_c0_join_fig.png", p2, width = 6.6, height = 6.2,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_join_fig.png")

## ---------- Figure 3: the crawl ---------------------------------------------
## Left panel: facet triangle (2-D) with gamma for decreasing delta.
Vt <- rbind(c(0, 0), c(1, 0), c(1/2, sqrt(3)/2))
projF <- function(G) { M <- G %*% Vt; data.frame(x = M[, 1], y = M[, 2]) }
sgc <- tan(seq(-1, 1, length.out = 3001) * (pi/2 - 0.004)) * 40
crawl_curve <- function(delta) {
  d <- c(2, -1 + delta, -1 - delta)
  data.frame(projF(softmax_rows(outer(sgc, d))), delta = delta)
}
cc <- do.call(rbind, lapply(c(0.4, 0.1, 0.02), crawl_curve))
cc$delta <- factor(cc$delta, levels = c(0.4, 0.1, 0.02))
## limit objects: symmetric arc (delta = 0) and the half-edge
arc <- data.frame(projF(softmax_rows(outer(sgc, c(2, -1, -1)))))
m23 <- (Vt[2, ] + Vt[3, ]) / 2
tri <- data.frame(x = Vt[c(1, 2, 3, 1), 1], y = Vt[c(1, 2, 3, 1), 2])
p3a <- ggplot() +
  geom_path(data = tri, aes(x, y), colour = "grey55", linewidth = 0.4) +
  geom_segment(aes(x = m23[1], y = m23[2], xend = Vt[3, 1], yend = Vt[3, 2]),
               colour = ORANGE, linewidth = 1.6, alpha = 0.9) +
  geom_path(data = arc, aes(x, y), colour = ORANGE, linewidth = 1.6,
            alpha = 0.9) +
  geom_path(data = cc, aes(x, y, group = delta, colour = delta),
            linewidth = 0.55) +
  scale_colour_manual(values = c("#BFD3D8", "#7FA3AC", TEAL),
                      name = expression(delta),
                      labels = c("0.4", "0.1", "0.02")) +
  geom_point(aes(x = m23[1], y = m23[2]), shape = 4, size = 2.6,
             stroke = 1.1) +
  annotate("text", x = Vt[1, 1] - 0.05, y = Vt[1, 2] - 0.03,
           label = "e[1]", parse = TRUE, size = 3.4) +
  annotate("text", x = Vt[2, 1] + 0.05, y = Vt[2, 2] - 0.03,
           label = "e[2]", parse = TRUE, size = 3.4) +
  annotate("text", x = Vt[3, 1], y = Vt[3, 2] + 0.045,
           label = "e[3]", parse = TRUE, size = 3.4) +
  annotate("text", x = m23[1] + 0.2, y = m23[2] + 0.06,
           label = "frac(1, 2)*(e[2] + e[3])", parse = TRUE, size = 3.1) +
  labs(title = "the log-line crawls as the tie break shrinks") +
  coord_equal(xlim = c(-0.08, 1.08), ylim = c(-0.06, 0.95)) +
  theme_void() +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        legend.position = c(0.86, 0.72))
## Right panel: limit surface in the tetrahedron = join over the arc
## (apex half) plus the gained half-face conv(e4, mid23, e3).
Gam0 <- gamma_curve(d0, sg)
Z0 <- rowSums(exp(outer(sg, d0))); t00 <- 1 / (1 + Z0)
cloud0 <- do.call(rbind, lapply(tt, function(t) {
  keep <- t >= t00
  if (!any(keep)) return(NULL)
  P <- cbind((1 - t) * Gam0[keep, 1:3, drop = FALSE], t)
  data.frame(proj4(P))
}))
gain <- rbind(E4[4, ], c(0, 1/2, 1/2, 0), c(0, 0, 1, 0))
gain2 <- data.frame(proj4(gain))
p3b <- ggplot() +
  geom_segment(data = tet_edges(), aes(x, y, xend = xend, yend = yend),
               colour = "grey75", linewidth = 0.35) +
  geom_polygon(data = gain2, aes(x, y), fill = ORANGE, alpha = 0.4,
               colour = ORANGE, linewidth = 0.5) +
  geom_point(data = cloud0, aes(x, y), colour = CLOUD, size = 0.3,
             alpha = 0.5) +
  geom_path(data = data.frame(proj4(Gam0)), aes(x, y), colour = TEAL,
            linewidth = 1.0) +
  geom_text(data = vlab, aes(x, y, label = lab), parse = TRUE, size = 3.8) +
  annotate("text", x = 0.72, y = -0.55,
           label = "gained~half*'-'*face:~sqrt(3)/4", parse = TRUE,
           size = 3.4, colour = ORANGE) +
  annotate("text", x = -0.62, y = -0.62,
           label = "apex~half~of~the~symmetric~join:~A[4]", parse = TRUE,
           size = 3.4, colour = "#5F8790") +
  labs(title = "the limit surface gains half a 2-face") +
  coord_equal(xlim = c(-1.2, 1.35)) + theme_void() +
  theme(plot.title = element_text(size = 11, hjust = 0.5))
p3 <- plot_grid(p3a, p3b, nrow = 1, labels = c("A", "B"), label_size = 12)
ggsave("../images/note_c0_crawl_fig.png", p3, width = 11.4, height = 5.4,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_crawl_fig.png")
