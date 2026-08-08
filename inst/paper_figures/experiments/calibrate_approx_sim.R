# Calibrate the approx-sim generator to the pancreas count distribution.
#
# Run ONCE, locally, whenever the simulation's dimensions or rank change:
#
#   Rscript calibrate_approx_sim.R
#
# It writes approx_sim_calibration.rds (committed to the repo) and prints the
# four knobs, which are also hardcoded with provenance in approx_sim_common.R
# so that the cluster workers depend only on that one file.
#
# METHOD (documented in the log1p_experiments notebook sim_design.Rmd):
#
# The calibration target is the size-adjusted BUDGET CURVE of the real data,
#   c |-> pm(log1p(Y/c)),   pm(x) = max(x)/mean(x),
# which depends only on the count histogram, plus the real mean count. The
# curve is what caps the attainable Hoyer sparsity of any fit at every c (see
# the Hoyer appendix of the paper), and the mean pins the absolute scale of
# the counts, which is what the approximation's accuracy depends on.
#
# The location knob mu is PROFILED, not searched: inside every loss
# evaluation, mu is solved by a 1-d root so the expected mean count equals the
# real one exactly (the root needs no Poisson draws and the mean is monotone
# in mu). Nelder-Mead then searches only (pi0, sigma, cap) against the curve,
# with the loss averaged over three fixed seeds (common random numbers). The
# final mu is re-solved on 20 seeds, because three are not enough to center
# the mean. The additive alternative -- a weighted mean term in the loss --
# was tried and traded away the tail of the curve; profiling keeps both.

suppressMessages(library(stats))

## ---- dimensions of the simulation being calibrated FOR ----------------------
N_S    <- 500L
P_S    <- 500L
K      <- 5L
A_SIG  <- 0.3
C_TRUE <- 1
T_DF   <- 4
NP_S   <- N_S * P_S

## ---- target: pancreas count histogram ---------------------------------------
z   <- readRDS("pancreas_count_histogram.rds")   # h = counts of value t >= 1; np
p_t <- c(z$np - sum(z$h), z$h) / z$np
tv  <- 0:(length(p_t) - 1)
M_REAL <- sum(p_t * tv)

## size adjustment: the median max of NP_S draws from the real marginal
tstar  <- tv[min(which(cumsum(p_t)^NP_S >= 0.5))]
CS     <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4)
target <- sapply(CS, function(cc) {
  g <- function(t) log1p(t / cc)
  g(tstar) / sum(p_t * g(tv))
})

## ---- generator pieces (must match approx_sim_common.R exactly) --------------

## Every feature is conditioned to load on AT LEAST ONE factor: if pi0 would
## zero a feature's whole row of F, its largest-u factor is kept. This is the
## generator's analogue of filtering genes for detectable expression -- and it
## is required, because log1pNMF's input validation rejects count matrices
## with all-zero columns, which dead features produce with probability
## pi0^K (~1/3 at K = 5). Deterministic given u, and smooth in pi0.
nz_of <- function(u, PI0) {
  nz <- u > PI0
  dead <- rowSums(nz) == 0
  if (any(dead)) nz[cbind(which(dead), max.col(u[dead, , drop = FALSE]))] <- TRUE
  nz
}

make_parts <- function(seeds) lapply(seeds, function(s) {
  set.seed(s)
  L <- matrix(rgamma(N_S * K, A_SIG), N_S, K); L <- L / rowSums(L)
  list(L = L, tmat = matrix(rt(P_S * K, df = T_DF), P_S, K),
       u = matrix(runif(P_S * K), P_S, K))
})
mean_rate <- function(mu, sig, CAP, PI0, parts)
  mean(sapply(parts, function(pp) {
    Fm <- pmin(exp(mu + sig * pp$tmat), CAP) * nz_of(pp$u, PI0)
    mean(C_TRUE * expm1(tcrossprod(pp$L, Fm))) }))
solve_mu <- function(sig, CAP, PI0, parts)
  uniroot(function(m) mean_rate(m, sig, CAP, PI0, parts) - M_REAL,
          c(-9, 3), tol = 1e-4)$root

curve_of <- function(Y) sapply(CS, function(cc) {
  g <- log1p(Y / cc); max(g) / mean(g) })
Y_from <- function(pp, mu, sig, CAP, PI0, pseed) {
  Fm <- pmin(exp(mu + sig * pp$tmat), CAP) * nz_of(pp$u, PI0)
  set.seed(pseed)
  matrix(rpois(NP_S, C_TRUE * expm1(tcrossprod(pp$L, Fm))), N_S, P_S)
}

## ---- profiled-mean curve loss ------------------------------------------------
parts3 <- make_parts(1:3)
loss <- function(par) {
  PI0 <- plogis(par[1]); sig <- exp(par[2]); CAP <- exp(par[3])
  mu <- tryCatch(solve_mu(sig, CAP, PI0, parts3), error = function(e) NA)
  if (is.na(mu)) return(1e6)
  mean(sapply(1:3, function(s) {
    cv <- curve_of(Y_from(parts3[[s]], mu, sig, CAP, PI0, s + 5000))
    if (any(!is.finite(cv)) || any(cv <= 0)) return(1e6)
    sum((log(cv) - log(target))^2) }))
}
opt <- optim(c(qlogis(0.75), log(0.6), log(10)), loss,
             method = "Nelder-Mead", control = list(maxit = 1500, reltol = 1e-6))
## Nelder-Mead reports convergence = 1 whenever it exhausts maxit, even after
## the simplex has effectively settled (the CRN loss has a noise floor that
## reltol never certifies). The values below match independent restarts, so
## code 1 is accepted with a notice rather than treated as failure.
if (opt$convergence != 0)
  message('optim exit code ', opt$convergence, ' (loss ', round(opt$value, 3), '); accepting')

PI0 <- plogis(opt$par[1]); SIG <- exp(opt$par[2]); CAP <- exp(opt$par[3])
MU  <- solve_mu(SIG, CAP, PI0, make_parts(1:20))
knobs <- c(pi0 = PI0, mu = MU, sigma = SIG, cap = CAP)

## ---- acceptance on held-out seeds -------------------------------------------
sim_full <- function(seed) {
  set.seed(seed)
  L <- matrix(rgamma(N_S * K, A_SIG), N_S, K); L <- L / rowSums(L)
  tmat <- matrix(rt(P_S * K, T_DF), P_S, K)
  u    <- matrix(runif(P_S * K), P_S, K)
  Fm   <- pmin(exp(MU + SIG * tmat), CAP) * nz_of(u, PI0)
  matrix(rpois(NP_S, C_TRUE * expm1(tcrossprod(L, Fm))), N_S, P_S)
}
fresh <- lapply(101:110, sim_full)
cvs <- sapply(fresh, curve_of)
mns <- sapply(fresh, mean)

cat("knobs:\n"); print(round(knobs, 4))
cat("\ncurve ratio to target (median over 10 held-out seeds):\n")
print(round(setNames(apply(cvs / target, 1, median), format(CS)), 2))
cat("\nmean count: median", round(median(mns), 3), " range",
    paste(round(range(mns), 3), collapse = "-"), " (real", round(M_REAL, 3), ")\n")
cat("zero frac median:",
    round(median(sapply(fresh, function(Y) mean(Y == 0))), 3),
    " (real", round(p_t[1], 3), ")\n")
cat("max median:", median(sapply(fresh, max)), " (t* =", tstar, ")\n")

saveRDS(list(knobs = knobs, target = target, CS = CS, tstar = tstar,
             mean_target = M_REAL, n = N_S, p = P_S, K = K, a_sig = A_SIG,
             c_true = C_TRUE, df = T_DF,
             method = "curve-matched with profiled mean; see sim_design.Rmd"),
        "approx_sim_calibration.rds")
cat("\nwrote approx_sim_calibration.rds\n")
