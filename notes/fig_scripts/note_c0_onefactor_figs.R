# Figures for the note's "one-factor story" section: log1p with an
# intercept plus ONE structural factor, vs additive models with at
# most two topics, at general p. Achievable sets are curves/segments;
# size = length.
#
#   Figure 1 (note_c0_staircase_fig.png): the cascade -- the general
#     construction, drawn at p = 4 where the simplex is a solid
#     tetrahedron. The curve for hierarchically tie-broken g =
#     (1.0101, 1.01, 1, 0) marches down the barycenter staircase
#     b4 (uniform) -> b3 -> b2 (edge midpoint) -> b1 = e1, with step
#     lengths 1/sqrt(j(j-1)).
#
#   Figure 2 (note_c0_lengths_fig.png): staircase length L(p) =
#     sum_{j=2}^p 1/sqrt(j(j-1)) against the two additive ceilings:
#     sqrt(2) (free K = 2: an edge) and sqrt(1 - 1/p) (rank-fair:
#     pinned at the uniform point). Crossover at p = 5; L(p) grows
#     like log p.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 11))

TEAL <- "#0E5C6B"; CLOUD <- "#9BB7BD"; ORANGE <- "#C2410C"

softmax_rows <- function(E) {
  E <- E - apply(E, 1, max)
  W <- exp(E); W / rowSums(W)
}

## ---------- Figure 1: the staircase (drawn at p = 4) ------------------------
T3 <- rbind(c( 1,  1, -1),
            c(-1, -1, -1),
            c( 1, -1,  1),
            c(-1,  1,  1)) / sqrt(3)
alignz <- function(v) {
  v <- v / sqrt(sum(v^2)); z <- c(0, 0, 1)
  ax <- c(v[2]*z[3]-v[3]*z[2], v[3]*z[1]-v[1]*z[3], v[1]*z[2]-v[2]*z[1])
  s <- sqrt(sum(ax^2)); cthet <- sum(v*z)
  if (s < 1e-12) return(diag(3))
  ax <- ax / s; K <- rbind(c(0,-ax[3],ax[2]), c(ax[3],0,-ax[1]),
                           c(-ax[2],ax[1],0))
  diag(3) + sqrt(1-cthet^2)*K + (1-cthet)*(K %*% K)
}
T3 <- T3 %*% t(alignz(T3[4, ] - colMeans(T3)))
az <- 100 * pi/180
Rz <- rbind(c(cos(az), -sin(az), 0), c(sin(az), cos(az), 0), c(0, 0, 1))
T3 <- T3 %*% t(Rz)
el <- 14 * pi/180
v1 <- c(1, 0, 0); v2 <- c(0, sin(el), cos(el))
proj4 <- function(L) {
  X <- L %*% T3
  data.frame(x = as.numeric(X %*% v1), y = as.numeric(X %*% v2))
}
E4 <- diag(4)
tet_edges <- do.call(rbind, apply(combn(4, 2), 2, function(k) {
  a <- proj4(E4[k[1], , drop = FALSE]); b <- proj4(E4[k[2], , drop = FALSE])
  data.frame(x = a$x, y = a$y, xend = b$x, yend = b$y)
}))
vlab4 <- data.frame(proj4(E4), lab = c("e[1]", "e[2]", "e[3]", "e[4]"))
vlab4$y <- vlab4$y + c(-0.09, -0.09, -0.09, 0.09)

g4 <- c(1.0101, 1.0100, 1.0000, 0)     # hierarchical gaps 1e-4, 1e-2, 1
s <- seq(0, 1, length.out = 60000)^2 * 6e5
P4 <- softmax_rows(outer(s, g4))
path4 <- proj4(P4)
bary <- rbind(rep(1/4, 4),
              c(1/3, 1/3, 1/3, 0),
              c(1/2, 1/2, 0, 0),
              c(1, 0, 0, 0))
bpts <- proj4(bary)
bpts$lab <- c("b[4] == uniform", "b[3]", "b[2]", NA)
bpts$dx <- c(0.28, 0.14, 0.11, 0)
bpts$dy <- c(0.10, 0.08, 0.09, 0)
seglab <- data.frame(
  x = (bpts$x[-4] + bpts$x[-1]) / 2 + c(-0.03, -0.06, -0.02),
  y = (bpts$y[-4] + bpts$y[-1]) / 2 + c(-0.10, -0.10, -0.10),
  lab = c("1/sqrt(12)", "1/sqrt(6)", "1/sqrt(2)"))
p1 <- ggplot() +
  geom_segment(data = tet_edges, aes(x, y, xend = xend, yend = yend),
               colour = "grey75", linewidth = 0.35) +
  geom_path(data = data.frame(x = bpts$x, y = bpts$y), aes(x, y),
            colour = CLOUD, linewidth = 1.6, alpha = 0.8,
            linetype = "22") +
  geom_path(data = path4, aes(x, y), colour = TEAL, linewidth = 1.0) +
  geom_point(data = bpts, aes(x, y), size = 2.2) +
  geom_text(data = bpts[1:3, ], aes(x + dx, y + dy, label = lab),
            parse = TRUE, size = 3.3) +
  annotate("text", x = bpts$x[4] - 0.04, y = bpts$y[4] - 0.16,
           label = "b[1] == e[1]", parse = TRUE, size = 3.3) +
  geom_text(data = seglab, aes(x, y, label = lab), parse = TRUE,
            size = 3.1, colour = "grey35") +
  geom_text(data = vlab4[2:4, ], aes(x, y, label = lab), parse = TRUE,
            size = 3.6) +
  labs(title = "the cascade, drawn at p = 4: down the barycenter staircase") +
  coord_equal(xlim = c(-1.2, 1.3)) + theme_void() +
  theme(plot.title = element_text(size = 11.5, hjust = 0.5))
ggsave("../images/note_c0_staircase_fig.png", p1, width = 7.6, height = 6.4,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_staircase_fig.png")

## ---------- Figure 2: lengths versus p --------------------------------------
pp <- 3:10
L <- sapply(pp, function(p) sum(1 / sqrt((2:p) * (1:(p-1)))))
df <- data.frame(p = pp, L = L)
p2 <- ggplot(df, aes(p, L)) +
  geom_hline(yintercept = sqrt(2), colour = "grey30", linewidth = 0.8) +
  geom_line(data = data.frame(p = seq(3, 10, 0.1)),
            aes(p, sqrt(1 - 1/p)), colour = "grey30", linewidth = 0.8,
            linetype = "22") +
  geom_line(colour = TEAL, linewidth = 0.8) +
  geom_point(colour = TEAL, size = 2.4) +
  annotate("text", x = 8.4, y = sqrt(2) - 0.06,
           label = "free~additive~ceiling*':'~an~edge*','~~sqrt(2)",
           parse = TRUE, size = 3.4, colour = "grey30") +
  annotate("text", x = 8.4, y = 1.03,
           label = "rank*'-'*fair~ceiling*':'~sqrt(1 - 1/p)",
           parse = TRUE, size = 3.4, colour = "grey30") +
  annotate("text", x = 4.7, y = 2.02,
           label = "staircase~length~L(p)", parse = TRUE, size = 3.6,
           colour = TEAL) +
  annotate("text", x = 5, y = 1.78, label = "crossover", size = 3.2,
           colour = ORANGE) +
  geom_point(aes(x = 5, y = L[3]), shape = 1, size = 4.5,
             colour = ORANGE, stroke = 1) +
  scale_x_continuous(breaks = 3:10) +
  labs(x = "p (number of features)", y = "length of achievable curve")
ggsave("../images/note_c0_lengths_fig.png", p2, width = 6.8, height = 4.2,
       dpi = 300, bg = "white")
message("Wrote ../images/note_c0_lengths_fig.png")
