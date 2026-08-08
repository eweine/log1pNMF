# Shared configuration and data generator for the exact-vs-approx simulation.
#
# Question: how much log-likelihood does the approximate objective give up
# relative to the exact one, on data with realistic count properties, across
# the c grid? Each dataset is generated ONCE per seed from the calibrated
# generator below, and every (c_fit, method) fits that same dataset from a
# rank-1 initialization until the objective converges.
#
# Sourced by run_approx_sim.R (the array worker).
#
# ---------------------------------------------------------------------------
# GENERATIVE MODEL (calibrated; see calibrate_approx_sim.R and the notebook
# sim_design.Rmd in the log1p_experiments repo for the full story)
#
#   L      n x K, rows ~ Dirichlet(A_SIG * 1_K): every sample concentrates on
#          a few factors; the loading sparsity is fixed by A_SIG and never
#          changes across seeds or c.
#   F      p x K, entries independently 0 with probability PI0, otherwise
#          exp(MU + SIGMA * t_4) truncated to [FLOOR, CAP]: a point mass at
#          zero and a heavy-tailed body, truncated on both sides. The CAP
#          bounds the maximum count (unbounded tails make max(Y)
#          seed-unstable; the cap is what the calibration aims at the real
#          data's size-adjusted maximum). The FLOOR is the counterpart of
#          expression filtering: without it, a feature active in a single
#          factor can draw a magnitude so small that its column's expected
#          total count is ~0.1 and Y gets an all-zero column no redraw can
#          fix (seed 7 did exactly this). At FLOOR = 0.05 the weakest
#          single-factor column has expected total count ~ 0.05 * n/K = 5.
#          The floor moves ~4% of the nonzero entries. Every feature is
#          additionally CONDITIONED to load on at least one factor, since
#          dead features (probability PI0^K, ~1/3 at K = 5) produce all-zero
#          columns deterministically, and log1pNMF rejects those.
#   theta  = L F'
#   lambda = C_TRUE * expm1(theta)
#   Y     ~ Poisson(lambda)         (no size factors; everything fits s=FALSE)
#
# The counts are redrawn (deterministically, seeded) in the rare event a row
# or column of Y is still all zero by Poisson chance.
#
# The four numbers (PI0, MU, SIGMA, CAP) come from calibrate_approx_sim.R,
# which matches the SIZE-ADJUSTED budget curve c |-> pm(log1p(Y/c)) of the
# real pancreas counts, with MU profiled so the mean count matches exactly.
#
# HONEST LIMITS of the calibration at n = p = 500, K = 5 (it was much tighter
# at K = 10: see sim_design.Rmd): with only 5 factors and Dirichlet-peaked
# rows, entry rates are nearly single-factor draws, and the knobs cannot feed
# the zero fraction, the mean, and the extreme tail all at once. Realized
# match on held-out seeds: mean count median 0.62 (real 0.674); zero fraction
# ~0.89 (real 0.822); budget curve within ~1.5x at the small-c end, 0.3-0.8x
# at the large-c end (median max ~2100 vs size-adjusted target 8820). Richer
# variants (free tail df, two-component magnitudes) fixed one end of the
# curve only by sacrificing the other, and were rejected.

## ---- configuration ---------------------------------------------------------
N_ROWS  <- 500L
N_COLS  <- 500L
K_TRUE  <- 5L
K_FIT   <- 5L

## the c grid the previous approximation-quality experiments used
CC_GRID <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4)
N_SEEDS <- 10L

## methods: the exact objective and both approximation techniques used in the
## paper's approximation-quality figure
METHODS <- c("exact", "taylor", "chebyshev")

## calibrated generator knobs -- output of calibrate_approx_sim.R
## (approx_sim_calibration.rds holds the same values plus the target curve).
## Provenance: pancreas counts (7606 x 18195), curve-matched with profiled
## mean at n = p = 500, K = 5, C_TRUE = 1, with the >= 1 factor conditioning.
A_SIG  <- 0.3
PI0    <- 0.7738
MU     <- -1.2512
SIGMA  <- 0.7067
CAP    <- 7.7587
FLOOR  <- 0.05    # fixed by design, not calibrated (see header)
T_DF   <- 4
C_TRUE <- 1

## fitting: run to convergence at TOL = 1e-6 (absolute change in the fit's own
## objective between successive iterations), rank-1 initialization always.
## MAXITER is a safety cap only (250k: far above anything the c-grid round
## needed), and there is deliberately NO optimization-time cap -- the SLURM
## wall clock (24h; see run_approx_sim.sbatch) is the sole backstop. See
## c_grid_sim_common.R for why a time cap is worse than none: it leaves fits
## stopped mid-descent, indistinguishable from a plateau.
MAXITER <- 250000L
TOL     <- 1e-6
THREADS <- 1L

OUT_DIR <- Sys.getenv("APPROX_SIM_OUT",
                      "/rafalab/eweine/log1p_experiments/approx_sim")

## ---- generative model ------------------------------------------------------

#' Nonzero pattern for F: entries survive with probability 1 - PI0, and any
#' feature that would end up with no factor keeps its largest-u one.
nz_of <- function(u, pi0 = PI0) {
  nz <- u > pi0
  dead <- rowSums(nz) == 0
  if (any(dead)) nz[cbind(which(dead), max.col(u[dead, , drop = FALSE]))] <- TRUE
  nz
}

#' The truth and counts for one seed. Deterministic given the seed, so every
#' array task sharing a seed rebuilds byte-identical data. The count draw is
#' retried (up to 50 sub-seeds) until no row or column of Y is entirely zero;
#' with the >= 1 factor conditioning this virtually always succeeds first try.
sim_dataset <- function(seed, n = N_ROWS, p = N_COLS, K = K_TRUE) {
  set.seed(seed)
  L  <- matrix(stats::rgamma(n * K, shape = A_SIG), n, K)
  L  <- L / rowSums(L)
  u  <- matrix(stats::runif(p * K), p, K)
  FF <- pmin(pmax(exp(MU + SIGMA * stats::rt(p * K, df = T_DF)), FLOOR), CAP) * nz_of(u)
  lam <- C_TRUE * expm1(tcrossprod(L, FF))

  for (try in 1:50) {
    set.seed(seed * 1000L + try)
    Y <- matrix(stats::rpois(n * p, lam), n, p)
    if (all(rowSums(Y) > 0) && all(colSums(Y) > 0)) break
  }
  if (any(rowSums(Y) == 0) || any(colSums(Y) == 0))
    stop("all-zero row/column after 50 redraws; seed ", seed)

  colnames(L) <- colnames(FF) <- paste0("k", seq_len(K))
  list(seed = seed, L = L, FF = FF, lambda = lam, Y = Y, y_draw = try)
}

pois_loglik <- function(Y, lam) {
  lam <- pmax(lam, 1e-300)
  sum(Y * log(lam) - lam - lgamma(Y + 1))
}

## ---- task id <-> (seed, cc, method) ----------------------------------------
## 0-based SLURM_ARRAY_TASK_ID over seeds x cc x method, method fastest. One
## task = ONE fit, so every fit is checkpointed independently and a hung node
## costs exactly one fit. 10 x 9 x 3 = 270 tasks, well under MaxArraySize.
N_TASKS <- N_SEEDS * length(CC_GRID) * length(METHODS)

task_spec <- function(id) {
  nm <- length(METHODS); nc <- length(CC_GRID)
  stopifnot(id >= 0, id < N_TASKS)
  meth_i <- id %% nm;  rest <- id %/% nm
  cc_i   <- rest %% nc
  seed_i <- rest %/% nc
  list(seed   = seed_i + 1L,
       cc     = CC_GRID[cc_i + 1L],
       method = METHODS[meth_i + 1L])
}

tag_of <- function(spec)
  sprintf("approxsim_seed%02d_c%s_%s", spec$seed,
          format(spec$cc, scientific = FALSE, trim = TRUE, drop0trailing = TRUE),
          spec$method)
