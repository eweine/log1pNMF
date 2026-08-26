# Shared configuration and data generator for the exact-vs-approx simulation,
# swept over the c used to GENERATE the data.
#
# Two questions share these runs. (1) Approximation quality: how much
# log-likelihood does the approximate objective give up relative to the exact
# one, across the fitting-c grid, on data with realistic count properties?
# (2) The simulation section: how do sparsity and related properties of the
# fits change with the fitting c, when the truth was generated at different
# c_true? Each dataset is generated once per (seed, c_true) from the
# calibrated generator below, and every (c_fit, method) fits that same
# dataset from a rank-1 initialization until the objective converges.
#
# Sourced by run_approx_sim.R (the array worker) and check_approx_sim.R.
#
# ---------------------------------------------------------------------------
# GENERATIVE MODEL (calibrated per c_true; see calibrate_approx_sim.R and the
# notebook sim_design.Rmd in the log1p_experiments repo)
#
#   L      n x K, rows ~ Dirichlet(A_SIG * 1_K): every sample concentrates on
#          a few factors. A_SIG never changes, so the true loading sparsity
#          is identical across seeds AND across c_true -- and because L
#          depends only on (seed, A_SIG), the same seed gives the SAME L at
#          every c_true.
#   F      p x K, entries independently 0 with probability pi0, otherwise
#          exp(mu + sigma * t_4) truncated to [floor, cap]. The cap bounds
#          the maximum count (unbounded tails make max(Y) seed-unstable);
#          the floor is the counterpart of expression filtering -- without
#          it, a feature active in a single factor can draw a magnitude so
#          small that its column of Y is all zero with probability ~1
#          (seed 7 did exactly this at c_true = 1), and log1pNMF rejects
#          all-zero columns. The floor is fixed on the RATE scale
#          (RATE_FLOOR expected counts) and mapped through each link, so it
#          means the same thing at every c_true. Every feature is
#          additionally CONDITIONED to load on at least one factor, since
#          dead features (probability pi0^K, ~1/3 at K = 5) produce all-zero
#          columns deterministically.
#          NOTE: balancing the F columns so w_k = mean(l_k)mean(f_k)/
#          mean(LF') = 1/K was implemented and then ABANDONED (see git
#          history around f4f1492): at c_true = Inf the identity link pins
#          every balanced column's mean to ~mean(rates)/(K mean(l_k)), so
#          max(Y) is structurally capped near ~300 vs the 8820 tail target
#          -- concentrated factor mass is how the identity link produces a
#          heavy tail. Rather than balance some regimes and not others,
#          the generator stays exactly the round-one law, so the round-one
#          outputs (seeds 1-10) remain valid and are NOT refit.
#   theta  = L F'
#   lambda = c_true * expm1(theta)   for finite c_true
#          = theta                   for c_true = Inf (plain Poisson NMF)
#   Y     ~ Poisson(lambda)          (no size factors; everything fits s=FALSE)
#
# The counts are redrawn (deterministically, seeded) in the rare event a row
# or column of Y is still all zero by Poisson chance.
#
# The knobs (pi0, mu, sigma, cap) are calibrated PER c_true against the SAME
# target -- the size-adjusted budget curve c |-> pm(log1p(Y/c)) of the real
# pancreas counts, plus the real mean count (profiled, matched exactly in
# expectation). One shared target is the principled version of "the datasets
# are comparable across c_true". See calibrate_approx_sim.R.
#
# HONEST LIMITS at n = p = 500, K = 5 (much tighter at K = 10: see
# sim_design.Rmd): with five factors and Dirichlet-peaked rows the knobs
# cannot feed the zero fraction, the mean, and the extreme tail all at once.
# Realized match: mean centered in expectation; zero fraction ~0.89 vs real
# 0.822; budget curve ~1.4-1.5x at the small-c end, ~0.3-0.8x at the large-c
# end. Richer variants (free tail df, two-component magnitudes) fixed one
# end only by sacrificing the other, and were rejected.

## ---- configuration ---------------------------------------------------------
N_ROWS  <- 500L
N_COLS  <- 500L
K_TRUE  <- 5L
K_FIT   <- 5L

## the sweep over the c used to generate the data
C_TRUE_GRID <- c(1e-3, 1, Inf)

## the c grid for FITTING: the round-one grid, PLUS c = Inf (plain Poisson
## NMF, i.e. the topic model, via fastTopics), fit with the EXACT objective
## only -- the quadratic zero-term approximations are constructions for the
## log1p objective and have no c = Inf analogue. Keeping the finite grid
## and the generator identical to round one means the 810 round-one fits
## (seeds 1-10) are reused as-is; only the 20 new seeds and the Inf fits
## run.
CC_GRID <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4, Inf)
N_SEEDS <- 30L

## The c_true = 1e-3 regime runs EXTRA seeds (31..N_SEEDS_SMALLC): its
## heavy-tailed local-optimum variation dominates the standard errors of
## the approximation-quality comparison, so it gets 100 seeds where the
## other regimes keep 30. The extra tasks are APPENDED to the spec table
## in seed order (31-60 first, then 61-100), so every previously assigned
## task id (and all existing outputs) is unchanged.
N_SEEDS_SMALLC <- 100L

## methods: the exact objective and both approximation techniques used in the
## paper's approximation-quality figure (finite c only; Inf is exact-only)
METHODS <- c("exact", "taylor", "chebyshev")

## calibrated generator knobs, one row per c_true -- output of
## calibrate_approx_sim.R (approx_sim_calibration.rds holds the same values
## plus the shared target curve). Provenance: pancreas counts (7606 x 18195),
## curve-matched with profiled mean at n = p = 500, K = 5. The floor is
## RATE_FLOOR = 0.05 expected counts mapped through each link.
A_SIG <- 0.3
T_DF  <- 4

## The knobs are read from the calibration file in this script's directory
## (the worker cd's there), so they cannot drift from what
## calibrate_approx_sim.R produced. Fail loudly if it is missing or stale.
cal_path <- "approx_sim_calibration.rds"
stopifnot(file.exists(cal_path))
CAL <- readRDS(cal_path)
stopifnot(identical(format(CAL$c_true_grid), format(C_TRUE_GRID)))
KNOBS <- do.call(rbind, lapply(CAL$per_c_true, function(r)
  data.frame(c_true = r$c_true, t(r$knobs))))
rownames(KNOBS) <- NULL

knobs_for <- function(cc) {
  i <- which(sapply(KNOBS$c_true, function(x) identical(x, cc) ||
                     (is.infinite(x) && is.infinite(cc)) ||
                     (is.finite(x) && is.finite(cc) && abs(x - cc) < 1e-12)))
  stopifnot(length(i) == 1)
  KNOBS[i, ]
}

## fitting: run to convergence at TOL = 1e-6 (absolute change in the fit's own
## objective between successive iterations), rank-1 initialization always.
## MAXITER is a safety cap only (1M -- round one showed taylor at small c
## needs up to ~850k), and there is deliberately NO optimization-time cap --
## the SLURM wall clock (24h; see run_approx_sim.sbatch) is the sole
## backstop. See c_grid_sim_common.R for why a time cap is worse than none.
## The c = Inf fits (fastTopics SCD) are capped at MAXITER_INF instead: SCD
## converges in thousands of iterations and its progress frame grows with
## the cap, so 1M would be pure waste there.
MAXITER     <- 1000000L
MAXITER_INF <- 100000L
TOL     <- 1e-6
THREADS <- 1L

## SAME directory as round one, on purpose: the generator and finite grid
## are unchanged, so the 810 existing outputs are this experiment's seeds
## 1-10 and the worker skips any task whose output already exists.
OUT_DIR <- Sys.getenv("APPROX_SIM_OUT",
                      "/rafalab/eweine/log1p_experiments/approx_sim")

## ---- generative model ------------------------------------------------------

#' Nonzero pattern for F: entries survive with probability 1 - pi0, and any
#' feature that would end up with no factor keeps its largest-u one.
nz_of <- function(u, pi0) {
  nz <- u > pi0
  dead <- rowSums(nz) == 0
  if (any(dead)) nz[cbind(which(dead), max.col(u[dead, , drop = FALSE]))] <- TRUE
  nz
}

#' The truth and counts for one (seed, c_true). Deterministic, so every array
#' task sharing (seed, c_true) rebuilds byte-identical data. L and the raw
#' draws behind F depend only on the seed, so the SAME L (and the same
#' uniforms/t-draws behind F) underlie every c_true; only the calibrated
#' magnitude law changes. The count draw is retried (up to 50 sub-seeds)
#' until no row or column of Y is entirely zero.
sim_dataset <- function(seed, c_true, n = N_ROWS, p = N_COLS, K = K_TRUE) {
  kb <- knobs_for(c_true)

  set.seed(seed)
  L  <- matrix(stats::rgamma(n * K, shape = A_SIG), n, K)
  L  <- L / rowSums(L)
  tm <- matrix(stats::rt(p * K, df = T_DF), p, K)
  u  <- matrix(stats::runif(p * K), p, K)
  FF <- pmin(pmax(exp(kb$mu + kb$sigma * tm), kb$floor), kb$cap) *
        nz_of(u, kb$pi0)

  th  <- tcrossprod(L, FF)
  lam <- if (is.infinite(c_true)) th else c_true * expm1(th)

  ci <- which(sapply(C_TRUE_GRID, function(x)
    (is.infinite(x) && is.infinite(c_true)) ||
    (is.finite(x) && is.finite(c_true) && abs(x - c_true) < 1e-12)))
  for (try in 1:50) {
    set.seed(seed * 1000L + 100L * ci + try)
    Y <- matrix(stats::rpois(n * p, lam), n, p)
    if (all(rowSums(Y) > 0) && all(colSums(Y) > 0)) break
  }
  if (any(rowSums(Y) == 0) || any(colSums(Y) == 0))
    stop("all-zero row/column after 50 redraws; seed ", seed,
         ", c_true ", format(c_true))

  colnames(L) <- colnames(FF) <- paste0("k", seq_len(K))
  list(seed = seed, c_true = c_true, L = L, FF = FF, lambda = lam, Y = Y,
       y_draw = try)
}

pois_loglik <- function(Y, lam) {
  lam <- pmax(lam, 1e-300)
  sum(Y * log(lam) - lam - lgamma(Y + 1))
}

## ---- task id <-> (seed, c_true, cc, method) --------------------------------
## The grid is no longer a full product (c = Inf is exact-only), so tasks are
## rows of an explicit spec table: method fastest, then cc, then c_true, then
## seed. One task = ONE fit, so every fit is checkpointed independently.
## 30 seeds x 3 c_true x (9 finite cc x 3 methods + Inf x exact) = 2520
## tasks, of which the 810 round-one fits already exist (the worker skips
## them). 2520 EXCEEDS the Slurm MaxArraySize default of 1001, so the grid
## is submitted in three slices via APPROX_SIM_TASK_OFFSET (see
## run_approx_sim.sbatch).
SPEC <- do.call(rbind, lapply(seq_len(N_SEEDS), function(si)
  do.call(rbind, lapply(C_TRUE_GRID, function(ct)
    do.call(rbind, lapply(CC_GRID, function(cc)
      data.frame(seed = si, c_true = ct, cc = cc,
                 method = if (is.infinite(cc)) "exact" else METHODS,
                 stringsAsFactors = FALSE)))))))
## extra small-c seeds, appended AFTER the original 2520 ids
if (N_SEEDS_SMALLC > N_SEEDS)
  SPEC <- rbind(SPEC, do.call(rbind,
    lapply((N_SEEDS + 1L):N_SEEDS_SMALLC, function(si)
      do.call(rbind, lapply(CC_GRID, function(cc)
        data.frame(seed = si, c_true = 1e-3, cc = cc,
                   method = if (is.infinite(cc)) "exact" else METHODS,
                   stringsAsFactors = FALSE))))))
rownames(SPEC) <- NULL
N_TASKS <- nrow(SPEC)

task_spec <- function(id) {
  stopifnot(id >= 0, id < N_TASKS)
  as.list(SPEC[id + 1L, ])
}

fmtc <- function(x) if (is.infinite(x)) "Inf" else
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

tag_of <- function(spec)
  sprintf("approxsim_seed%02d_ctrue%s_cfit%s_%s", spec$seed,
          fmtc(spec$c_true), fmtc(spec$cc), spec$method)
