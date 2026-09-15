# Numbers quoted in the "area supremum for p >= 4" subsection of
# note_c0_geometry.tex (no figure output).
#
# Setting: p = 4, two structural factors plus an intercept. The
# achievable set is sigma(H) for a half-plane H of a 2-plane V in the
# centered space; area is 2-D surface (Hausdorff) measure inside the
# 3-simplex.
#
# (1) half_area: direct quadrature of the surface integral over a
#     half-plane, with the area element computed via Cauchy-Binet
#     (sum of squared 2x2 minors -- the naive Gram determinant
#     n1*n2 - cr^2 cancels catastrophically where the two tangent
#     vectors are near-parallel and must not be used).
#     Used for the "embedded p = 3" candidate: normal c = centered
#     e4, surface {lam4^3 = lam1 lam2 lam3}, half {lam3 >= sqrt(lam1
#     lam2)}: area 0.7567.
#
# (2) join_areas: exact 1-D formulas for the join surfaces of the
#     note's Lemma (normal supported on 3 coordinates): sigma(V) is
#     the cone with apex e4 over the log-line gamma(s) =
#     softmax(s*d) in the opposite facet, a ruled surface with area
#     element (1-t)*w(s).
#
# (3) sweep of the join direction d(theta) toward the tie at theta =
#     30 deg (d = (2,-1,-1)): the apex-half area converges to
#     0.4561 + sqrt(3)/4 = 0.8891 (limit surface gains the half
#     facet conv{e4, mid(e2,e3), e3}) while the value AT the tie is
#     0.4561 -- the discontinuity/non-attainment phenomenon.

## ---- (1) direct half-plane surface quadrature ------------------------------
half_area <- function(w, u, ns = 480, nt = 240) {
  p <- length(w)
  s <- (seq_len(ns) - 0.5) / ns * pi - pi/2      # a = tan(s)
  t <- (seq_len(nt) - 0.5) / nt * (pi/2)         # b = tan(t) >= 0
  g <- expand.grid(s = s, t = t)
  a <- tan(g$s); b <- tan(g$t)
  wt <- (1/cos(g$s)^2) * (1/cos(g$t)^2) * (pi/ns) * (pi/2/nt)
  eta <- outer(a, w) + outer(b, u)
  eta <- eta - apply(eta, 1, max)
  S <- exp(eta); S <- S / rowSums(S)
  j1 <- S * (matrix(w, nrow(S), p, byrow = TRUE) - c(S %*% w))
  j2 <- S * (matrix(u, nrow(S), p, byrow = TRUE) - c(S %*% u))
  d2 <- 0
  for (i in 1:(p-1)) for (j in (i+1):p)
    d2 <- d2 + (j1[, i] * j2[, j] - j1[, j] * j2[, i])^2
  sum(sqrt(d2) * wt)
}
onb <- function(w, u) {
  w <- w - mean(w); w <- w / sqrt(sum(w^2))
  u <- u - mean(u); u <- u - sum(u * w) * w
  list(w = w, u = u / sqrt(sum(u^2)))
}
f <- onb(c(1, -1, 0, 0), c(-1, -1, 2, 0))   # embedded p = 3 optimum
cat(sprintf("embedded p = 3 candidate: area = %.4f\n",
            half_area(f$w, f$u)))

## ---- (2) exact join areas (1-D integrals) ----------------------------------
join_areas <- function(d, n = 60000) {
  d <- d - mean(d)
  x <- (seq_len(n) - 0.5) / n * pi - pi/2
  s <- tan(x); jac <- 1/cos(x)^2 * (pi/n)
  eta <- outer(s, d); m <- apply(eta, 1, max)
  E <- exp(eta - m); Zr <- rowSums(E); G <- E / Zr
  logZ <- m + log(Zr)                            # log sum_j exp(s d_j)
  Gd <- c(G %*% d)
  Gp <- G * (matrix(d, n, 3, byrow = TRUE) - Gd) # Gamma'(s), facet part
  A <- cbind(Gp, 0); B <- cbind(-G, 1)           # e4 - Gamma = (-gamma, 1)
  w <- sqrt(pmax(rowSums(A^2) * rowSums(B^2) - rowSums(A*B)^2, 0))
  t0 <- plogis(-logZ)                            # 1/(1 + Z)
  c(full = sum(w * jac) / 2,
    facet_half = sum(w * (t0 - t0^2/2) * jac),
    apex_half  = sum(w * (1 - t0)^2/2 * jac))
}
for (d in list(c(1, 0, -1), c(2, -1, -1))) {
  a <- join_areas(d)
  cat(sprintf("join d = (%s): full %.4f | facet half %.4f | apex half %.4f\n",
              paste(d, collapse = ","), a[1], a[2], a[3]))
}

## ---- (3) discontinuity at the tie direction --------------------------------
b1 <- c(1, -1, 0) / sqrt(2); b2 <- c(1, 1, -2) / sqrt(6)
cat("\napex-half area approaching the tie at theta = 30 deg:\n")
for (dth in c(2, 1, 0.5, 0.2, 0.1, 0.05, 0.02)) {
  th <- (30 + dth) * pi / 180
  a <- join_areas(cos(th) * b1 + sin(th) * b2)
  cat(sprintf("  theta = 30 + %-5g deg: apex half = %.4f\n", dth, a["apex_half"]))
}
cat(sprintf("limit-complex value 0.4561 + sqrt(3)/4 = %.4f\n",
            0.4561 + sqrt(3)/4))
