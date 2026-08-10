# Calibrate the approx-sim generator to the pancreas count distribution, ONCE
# PER c_true. Run locally whenever the dimensions, rank, or c_true grid change:
#
#   Rscript calibrate_approx_sim.R
#
# It writes approx_sim_calibration.rds (a list keyed by c_true, committed to
# the repo) and prints the knobs, which are also hardcoded with provenance in
# approx_sim_common.R so the cluster workers depend only on that one file.
#
# METHOD (documented in the log1p_experiments notebook sim_design.Rmd):
#
# The calibration target -- IDENTICAL for every c_true, which is the whole
# point -- is the size-adjusted BUDGET CURVE of the real data,
#   c |-> pm(log1p(Y/c)),   pm(x) = max(x)/mean(x),
# plus the real mean count. The curve depends only on the count histogram and
# governs the factor sparsity attainable at each fitting c (see the Hoyer
# appendix of the paper); the mean pins the absolute scale of the counts,
# which is what the approximation's accuracy depends on. Solving the knobs
# per c_true against one shared target is the principled version of "the
# datasets are comparable across c_true".
#
# The location knob mu is PROFILED, not searched: inside every loss
# evaluation it is solved by a 1-d root so the expected mean count equals the
# real one exactly. Nelder-Mead searches only (pi0, sigma, cap) against the
# curve, with the loss averaged over three fixed seeds (common random
# numbers); the final mu is re-solved on 20 seeds. Nelder-Mead's exit code 1
# (maxit) is accepted: the CRN loss has a noise floor that reltol never
# certifies, and restarts reproduce the same knobs.

suppressMessages(library(stats))

## ---- dimensions of the simulation being calibrated FOR ----------------------
N_S    <- 500L
P_S    <- 500L
K      <- 5L
A_SIG  <- 0.3
T_DF   <- 4
NP_S   <- N_S * P_S

C_TRUE_GRID <- c(1e-3, 1, Inf)

## the magnitude floor is fixed on the RATE scale (a retained feature has
## detectable expression; see approx_sim_common.R for the dead-column
## rationale) and mapped through each link:
RATE_FLOOR <- 0.05
th_of  <- function(r, cc)  if (is.infinite(cc)) r  else log1p(r / cc)
lam_of <- function(th, cc) if (is.infinite(cc)) th else cc * expm1(th)

## ---- target: pancreas count histogram (shared by every c_true) --------------
z   <- readRDS("pancreas_count_histogram.rds")   # h = counts of value t >= 1; np
p_t <- c(z$np - sum(z$h), z$h) / z$np
tv  <- 0:(length(p_t) - 1)
M_REAL <- sum(p_t * tv)

tstar  <- tv[min(which(cumsum(p_t)^NP_S >= 0.5))]  # median max of NP_S draws
CS     <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4)
target <- sapply(CS, function(cc) {
  g <- function(t) log1p(t / cc)
  g(tstar) / sum(p_t * g(tv))
})

## ---- generator pieces (must match approx_sim_common.R exactly) --------------
nz_of <- function(u, PI0) {
  nz <- u > PI0
  dead <- rowSums(nz) == 0
  if (any(dead)) nz[cbind(which(dead), max.col(u[dead, , drop = FALSE]))] <- TRUE
  nz
}
BAL_TOL  <- 1e-3
BAL_ITER <- 200L
balance_F <- function(L, FF, nz, FLOOR, CAP) {
  KK <- ncol(FF)
  Lm <- colMeans(L)
  for (i in seq_len(BAL_ITER)) {
    w <- Lm * colMeans(FF); w <- w / sum(w)
    if (max(abs(KK * w - 1)) < BAL_TOL) break
    FF <- sweep(FF, 2, (1 / KK) / w, "*")
    FF <- pmin(pmax(FF, FLOOR), CAP) * nz
  }
  FF
}
make_parts <- function(seeds) lapply(seeds, function(s) {
  set.seed(s)
  L <- matrix(rgamma(N_S * K, A_SIG), N_S, K); L <- L / rowSums(L)
  list(L = L, tmat = matrix(rt(P_S * K, df = T_DF), P_S, K),
       u = matrix(runif(P_S * K), P_S, K))
})
## BAL: balance only at finite c_true. At the identity link, balance caps
## max(F) at p * mean(f_k) ~ 337 (tail unmatchable, curve loss ~16.5 from
## every start) and pushes the mass onto rare cap-sized draws (held-out
## mean counts 0.016-1.1) -- concentrated factor mass is how the identity
## link produces a heavy tail, so c_true = Inf keeps the unbalanced law.
Fm_of <- function(pp, mu, sig, CAP, PI0, FLOOR, BAL) {
  nz <- nz_of(pp$u, PI0)
  F0 <- pmin(pmax(exp(mu + sig * pp$tmat), FLOOR), CAP) * nz
  if (BAL) F0 <- balance_F(pp$L, F0, nz, FLOOR, CAP)
  F0
}
curve_of <- function(Y) sapply(CS, function(cc) {
  g <- log1p(Y / cc); max(g) / mean(g) })

calibrate_one <- function(cc) {
  FLOOR <- th_of(RATE_FLOOR, cc)
  BAL   <- is.finite(cc)
  mean_rate <- function(mu, sig, CAP, PI0, parts)
    mean(sapply(parts, function(pp)
      mean(lam_of(tcrossprod(pp$L, Fm_of(pp, mu, sig, CAP, PI0, FLOOR, BAL)), cc))))
  solve_mu <- function(sig, CAP, PI0, parts)
    uniroot(function(m) mean_rate(m, sig, CAP, PI0, parts) - M_REAL,
            c(log(FLOOR) - 10, log(CAP) + 0.5), tol = 1e-4)$root
  Y_from <- function(pp, mu, sig, CAP, PI0, pseed) {
    lam <- lam_of(tcrossprod(pp$L, Fm_of(pp, mu, sig, CAP, PI0, FLOOR, BAL)), cc)
    set.seed(pseed)
    matrix(rpois(NP_S, lam), N_S, P_S)
  }

  parts3 <- make_parts(1:3)
  loss <- function(par) {
    PI0 <- plogis(par[1]); sig <- exp(par[2]); CAP <- exp(par[3])
    if (CAP <= FLOOR) return(1e6)
    mu <- tryCatch(solve_mu(sig, CAP, PI0, parts3), error = function(e) NA)
    if (is.na(mu)) return(1e6)
    mean(sapply(1:3, function(s) {
      cv <- curve_of(Y_from(parts3[[s]], mu, sig, CAP, PI0, s + 5000))
      if (any(!is.finite(cv)) || any(cv <= 0)) return(1e6)
      sum((log(cv) - log(target))^2) }))
  }

  ## c_true-aware starting values via the inverse link: a typical nonzero
  ## entry should contribute a rate ~0.3, and the cap should reach the rate
  ## scale of the previous round's realized maximum (~2000).
  init <- c(qlogis(0.77), log(0.7), log(th_of(2000, cc)))
  init[3] <- max(init[3], log(FLOOR * 3))   # keep cap > floor at small c_true
  ## MULTI-START: balancing the factor columns reshaped the loss surface,
  ## and from the generic start alone Nelder-Mead degenerated at
  ## c_true = Inf (mu ~ -13, all mass at the clamps). The knobs of the
  ## UNBALANCED round-one calibration are a second start; keep the better
  ## solution.
  hints <- list("0.001" = c(0.9592, 0.0994, 15.5685),   # (pi0, sigma, cap)
                "1"     = c(0.8379, 0.6652, 7.8197),
                "Inf"   = c(0.9011, 1.4668, 911.2182))
  h <- hints[[format(cc)]]
  starts <- list(init, c(qlogis(h[1]), log(h[2]), log(h[3])))
  opts <- lapply(starts, function(s)
    optim(s, loss, method = "Nelder-Mead",
          control = list(maxit = 800, reltol = 1e-6)))
  vals <- sapply(opts, function(o) o$value)
  message("  start losses: ", paste(round(vals, 3), collapse = " / "))

  ## acceptance on held-out seeds, per candidate
  sim_full <- function(seed, MU, SIG, CAP, PI0) {
    set.seed(seed)
    L  <- matrix(rgamma(N_S * K, A_SIG), N_S, K); L <- L / rowSums(L)
    tm <- matrix(rt(P_S * K, T_DF), P_S, K)
    u  <- matrix(runif(P_S * K), P_S, K)
    nz <- nz_of(u, PI0)
    Fm <- pmin(pmax(exp(MU + SIG * tm), FLOOR), CAP) * nz
    if (BAL) Fm <- balance_F(L, Fm, nz, FLOOR, CAP)
    matrix(rpois(NP_S, lam_of(tcrossprod(L, Fm), cc)), N_S, P_S)
  }
  sols <- lapply(opts, function(o) {
    PI0 <- plogis(o$par[1]); SIG <- exp(o$par[2]); CAP <- exp(o$par[3])
    MU  <- solve_mu(SIG, CAP, PI0, make_parts(1:20))
    fresh <- lapply(101:110, sim_full, MU = MU, SIG = SIG, CAP = CAP, PI0 = PI0)
    mns <- sapply(fresh, mean)
    list(value = o$value, convergence = o$convergence,
         knobs = c(pi0 = PI0, mu = MU, sigma = SIG, cap = CAP, floor = FLOOR),
         fresh = fresh, mns = mns, spread = diff(range(mns)))
  })

  ## Pick by loss, BUT near-tied losses go to the stabler generator: the
  ## balanced c_true = Inf surface has near-tied optima whose realized mean
  ## counts differ 10x in seed-to-seed spread, and a 30-seed experiment
  ## cannot afford datasets whose total mass swings by a factor of 70.
  pick <- which.min(vals)
  for (j in seq_along(sols))
    if (j != pick && sols[[j]]$value <= 1.05 * sols[[pick]]$value &&
        sols[[j]]$spread < 0.5 * sols[[pick]]$spread) pick <- j
  message("  candidate mean-count spreads: ",
          paste(round(sapply(sols, `[[`, "spread"), 3), collapse = " / "),
          "; picked start ", pick)
  sol <- sols[[pick]]
  if (sol$convergence != 0)
    message("  optim exit ", sol$convergence, " (loss ",
            round(sol$value, 3), "); accepting")

  knobs <- sol$knobs
  fresh <- sol$fresh
  cvs <- sapply(fresh, curve_of)
  mns <- sol$mns
  opt <- list(value = sol$value)

  cat(sprintf("\n===== c_true = %s   (loss %.3f)\n", format(cc), opt$value))
  print(round(knobs, 4))
  cat("curve ratio to target (median over 10 held-out seeds):\n")
  print(round(setNames(apply(cvs / target, 1, median), format(CS)), 2))
  cat("mean count: median", round(median(mns), 3), " range",
      paste(round(range(mns), 3), collapse = "-"), " (real", round(M_REAL, 3), ")\n")
  cat("zero frac median:",
      round(median(sapply(fresh, function(Y) mean(Y == 0))), 3),
      " (real", round(p_t[1], 3), ")   max median:",
      median(sapply(fresh, max)), " (t* =", tstar, ")\n")

  list(c_true = cc, knobs = knobs, loss = opt$value)
}

res <- lapply(C_TRUE_GRID, calibrate_one)
names(res) <- format(C_TRUE_GRID)

saveRDS(list(per_c_true = res, target = target, CS = CS, tstar = tstar,
             mean_target = M_REAL, n = N_S, p = P_S, K = K, a_sig = A_SIG,
             df = T_DF, rate_floor = RATE_FLOOR,
             c_true_grid = C_TRUE_GRID,
             method = "curve-matched with profiled mean, per c_true; see sim_design.Rmd"),
        "approx_sim_calibration.rds")
cat("\nwrote approx_sim_calibration.rds\n")
