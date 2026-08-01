# Shared configuration, data generator and metrics for the c-grid simulation.
#
# Question: how does the log1p model behave when the *fitted* c differs from the
# c that generated the data? We simulate rank-K data at each c_true on a grid,
# and fit the model at every c_fit on the same grid -- a full cross of link
# parameters -- over several independent seeds.
#
# Sourced by run_c_grid_sim.R (the array worker) and by
# check_c_grid_sim_data.R (the local no-fitting diagnostic).
#
# ---------------------------------------------------------------------------
# GENERATIVE MODEL
#
#   L      n x K, rows ~ Dirichlet(a_sig * 1_K), a_sig < 1 so the rows are
#          peaked toward the pure states. Rows sum to 1 -- this is what makes
#          the offset trick below exact, and the near-pure rows are what make
#          the factorization identifiable.
#   F0     p x K, sparse: each entry is 0 with probability pi0, otherwise
#          LOG-UNIFORM over a range of F_RANGE. The log-uniform magnitudes are
#          what give features a wide dynamic range; there is no separate
#          baseline factor, so no column can dominate the theta budget.
#
#          Magnitudes are log-uniform (bounded) rather than lognormal on
#          purpose. At small c the link is exponential, so an unbounded right
#          tail in theta0 gets exponentiated: with lognormal magnitudes the
#          c = 1e-3 datasets came out with max(Y) ~ 1e16 while the c = Inf ones
#          had max(Y) ~ 90. Bounding them is what makes max(Y) comparable
#          across the grid at all.
#   theta0 = L F0'                                            (n x p)
#   theta  = alpha_c * theta0 + beta_c
#   lambda = c * expm1(theta)          for finite c
#          = theta                     for c = Inf   (the c -> Inf limit, i.e.
#                                                     plain Poisson NMF)
#   Y_ij  ~ Poisson(lambda_ij)          (no size factors; everything fits s=FALSE)
#
# Because rowSums(L) == 1, theta = alpha*theta0 + beta is exactly L %*% t(F_c)
# with F_c = alpha*F0 + beta. So the realized truth is still a rank-K
# nonnegative factorization -- the offset costs no rank. F_c is not itself
# sparse (it has a floor of beta); the *sparse* object is F0, and both are
# scored in the metrics below.
#
# ---------------------------------------------------------------------------
# WHY TWO KNOBS
#
# With only a scale on F, the datasets would not be remotely comparable across
# c. At small c the link is log-like, so a fixed spread in theta becomes a
# multiplicative spread in lambda and max(lambda) explodes; at large c the link
# is linear and the spread stays proportional. (alpha_c, beta_c) are therefore
# solved per c so that every dataset hits the same two targets:
#
#   target_zero  the fraction of zero entries in Y,  E[exp(-lambda)]
#   target_q     the qprob quantile of lambda       (controls dynamic range)
#
# The solve is two nested uniroots and is exactly monotone in each direction:
# zero-fraction decreases in beta at fixed alpha, and the upper quantile
# increases in alpha once beta has been re-solved. It is deterministic given
# the seed, so every array task reconstructs an identical dataset.

## ---- configuration ---------------------------------------------------------
N_ROWS      <- 300L
N_COLS      <- 300L
K_TRUE      <- 3L
K_FIT       <- 3L

C_GRID      <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, Inf)
N_SEEDS     <- 10L
INITS       <- c("rank1", "random")   # every cell is fit from both

## generative shape
A_SIG       <- 0.3     # Dirichlet concentration (< 1 => rows peaked / sparse)
PI0         <- 0.5     # P(F0[j, k] == 0)
F_RANGE     <- 1e3     # nonzero F0 magnitudes span this factor, log-uniform

## calibration targets (identical for every c, which is the whole point).
## Tuned so that all eight datasets land on the same zero fraction and the same
## upper tail; see check_c_grid_sim_data.R for the realized spread.
TARGET_ZERO <- 0.50    # fraction of zeros in Y
TARGET_Q    <- 20      # QPROB-quantile of lambda
QPROB       <- 0.999

## fitting: run to convergence, with maxiter and max_time only as safety caps.
##
## TOL is an ABSOLUTE log-likelihood difference between successive iterations.
## The log-likelihood here is ~2.8e5, so double precision bottoms out around
## 2.8e5 * 1e-16 = 3e-11: a tolerance of 1e-10 can essentially never fire, which
## is why every fit at that setting reported converged = FALSE even after the
## objective had stopped moving entirely. 1e-6 is still utterly negligible
## against 2.8e5 (a relative change of 4e-12) and does fire.
## Cost is very uneven across the grid: c_fit >= 1 converges in a few hundred
## iterations, while the small-c_fit corner needs thousands and the
## worst-misspecified corner (linear-scale truth fit with a near-log link) is the
## slowest of all. There is deliberately NO optimization-time cap now -- an
## earlier round used one and left cells stopped mid-descent, which is
## indistinguishable in the output from a cell that genuinely plateaued.
## MAXITER is the only stop besides TOL, and the SLURM wall clock is the backstop.
MAXITER     <- 100000L
TOL         <- 1e-6
MAX_TIME    <- Inf     # no optimization-time cap; MAXITER / TOL stop the fit
THREADS     <- 1L      # one core per array task

OUT_DIR     <- Sys.getenv("C_GRID_OUT",
                          "/rafalab/eweine/log1p_experiments/c_grid_sim")

## rerun bookkeeping (see plan_c_grid_rerun.R / run_c_grid_rerun.R)
WORKLIST    <- file.path(OUT_DIR, "c_grid_rerun_todo.rds")
SUPERSEDED  <- file.path(OUT_DIR, "superseded")

## ---- generative model ------------------------------------------------------

#' Draw the shared truth (L, F0) for one seed. Identical across all c.
sim_truth <- function(seed, n = N_ROWS, p = N_COLS, K = K_TRUE) {
  set.seed(seed)

  G <- matrix(stats::rgamma(n * K, shape = A_SIG), n, K)
  L <- G / rowSums(G)

  lr <- log(F_RANGE) / 2
  F0 <- matrix(0, p, K)
  for (k in seq_len(K)) {
    nz <- stats::runif(p) > PI0
    F0[nz, k] <- exp(stats::runif(sum(nz), -lr, lr))
  }

  colnames(L) <- colnames(F0) <- paste0("k", seq_len(K))
  rownames(L) <- paste0("i", seq_len(n))
  rownames(F0) <- paste0("j", seq_len(p))
  list(L = L, F0 = F0, theta0 = L %*% t(F0))
}

#' Rates under (alpha, beta) at a given c. c = Inf is the linear limit.
lambda_of <- function(theta0, alpha, beta, cc) {
  th <- alpha * theta0 + beta
  if (is.infinite(cc)) return(th)
  cc * expm1(pmin(th, 700))          # pmin guards the c -> 0 corner
}

zero_frac <- function(lam) mean(exp(-lam))

#' Solve (alpha_c, beta_c) so that lambda hits TARGET_ZERO and TARGET_Q.
#'
#' Two nested uniroots. Brackets are derived from a first-order guess rather
#' than fixed wide ones: a fixed bracket in log(alpha) is unusable across six
#' decades of c, because at the small-c end a large alpha saturates the link and
#' the inner problem stops having a root at all.
calibrate <- function(theta0, cc,
                      target_zero = TARGET_ZERO, target_q = TARGET_Q,
                      qprob = QPROB) {

  ## theta that produces a given rate, i.e. the inverse link
  th_of_lam <- function(l) if (is.infinite(cc)) l else log1p(l / cc)

  ## a constant rate lam_mid would by itself give the target zero fraction
  lam_mid <- -log(target_zero)
  th_hi   <- th_of_lam(target_q)
  th_mid  <- th_of_lam(lam_mid)
  beta_hi <- th_of_lam(1e6)            # ample upper bound: rate 1e6 is absurd

  ## inner: zero-fraction decreases in beta. beta >= 0 is required (F_c =
  ## alpha*F0 + beta must stay nonnegative), so if alpha is already so large
  ## that beta = 0 overshoots, clamp -- the outer search then moves alpha down.
  beta_for <- function(alpha) {
    f <- function(b) zero_frac(lambda_of(theta0, alpha, b, cc)) - target_zero
    if (f(0) <= 0) return(0)
    stats::uniroot(f, c(0, beta_hi), tol = 1e-12)$root
  }

  ## outer: the upper quantile increases in alpha once beta is re-solved.
  ## The bracket is capped at alpha * max(theta0) = 200 so the link can never
  ## saturate: beyond that every entry sits in the flat tail, the quantile stops
  ## being monotone in alpha, and the sign change disappears. Below the cap both
  ## ends are guaranteed: alpha -> 0 gives a near-constant rate (quantile below
  ## target), alpha at the cap gives an enormous one.
  g <- function(la) {
    alpha <- exp(la)
    lam   <- lambda_of(theta0, alpha, beta_for(alpha), cc)
    as.numeric(stats::quantile(lam, qprob)) - target_q
  }
  la_hi <- log(200 / max(theta0))
  alpha <- exp(stats::uniroot(g, c(la_hi - 30, la_hi), tol = 1e-10)$root)

  list(alpha = alpha, beta = beta_for(alpha))
}

#' Full dataset for one (seed, c_true): truth, calibration, rates, counts.
#' Deterministic, so every array task sharing (seed, c_true) rebuilds the same Y.
sim_dataset <- function(seed, cc, truth = NULL) {
  if (is.null(truth)) truth <- sim_truth(seed)
  cal <- calibrate(truth$theta0, cc)

  Fc  <- cal$alpha * truth$F0 + cal$beta
  lam <- lambda_of(truth$theta0, cal$alpha, cal$beta, cc)

  ## Y seed is a deterministic function of (seed, c_true) and independent of the
  ## truth draw, so all 8 c_fit tasks see byte-identical counts.
  set.seed(seed * 1000L + match(cc, C_GRID))
  Y <- matrix(stats::rpois(length(lam), lam), nrow(lam), ncol(lam))
  dimnames(Y) <- dimnames(lam) <- list(rownames(truth$L), rownames(truth$F0))

  list(truth = truth, cc = cc, seed = seed, alpha = cal$alpha, beta = cal$beta,
       FF_true = Fc, lambda = lam, Y = Y)
}

#' Marginal summaries used to check that datasets are comparable across c.
data_summary <- function(d) {
  data.frame(
    c_true    = d$cc,
    alpha     = d$alpha,
    beta      = d$beta,
    zero_frac = mean(d$Y == 0),
    mean_Y    = mean(d$Y),
    max_Y     = max(d$Y),
    q999_Y    = as.numeric(stats::quantile(d$Y, 0.999)),
    mean_lam  = mean(d$lambda),
    med_lam   = stats::median(d$lambda),
    q999_lam  = as.numeric(stats::quantile(d$lambda, 0.999)),
    max_lam   = max(d$lambda),
    # spread of per-column (gene) mean rate: the dynamic range across features
    gene_mean_q10 = as.numeric(stats::quantile(colMeans(d$lambda), 0.10)),
    gene_mean_q90 = as.numeric(stats::quantile(colMeans(d$lambda), 0.90)))
}

## ---- fitting ---------------------------------------------------------------

RANK1_PAD <- 1e-8   # value the rank-1 warm start pads the other K-1 columns with
                    # (matches what fit_poisson_log1p_nmf's own rank1 path uses)

## The random initialization, matching what fit_poisson_log1p_nmf draws
## internally for init_method = "random". Building it explicitly means the
## c_fit = Inf column can be given the SAME random start as every other column,
## instead of fastTopics' own scheme -- otherwise "random" would not mean the
## same thing across the grid.
random_init <- function(n, p, K) {
  list(LL = matrix(stats::runif(n * K, 1e-8, 0.05), n, K),
       FF = matrix(stats::runif(p * K, 1e-8, 0.05), p, K))
}

#' Fit one cell.
#'
#' `init` is "rank1" or "random", and means the same thing at every c_fit
#' including Inf. For rank1 at c_fit = Inf the analogue is fastTopics' exact
#' rank-1 Poisson NMF solution padded to K, following the same pattern as
#' factor_sim_identifiable_reduc_fig.R; for random it is the same runif draw the
#' log1p package uses. The first round of this experiment left the c_fit = Inf
#' column on fastTopics' default `topicscore` initialization, which made that
#' column not init-comparable to the rest of the grid.
#'
#' Warm start (init_LL / init_FF given) resumes an earlier fit and overrides
#' `init`. For the log1p path this is exact: the package divides the supplied
#' init by sqrt(max(1, cc)) on the way in and multiplies back on the way out, so
#' feeding back a previous fit's returned LL / FF at the SAME cc reproduces its
#' internal theta bit for bit. Verified: 60 iterations then a 340-iteration warm
#' continuation gives an identical log-likelihood to a single cold 400-iteration
#' run.
fit_cell <- function(Y, cc, init = "rank1", init_LL = NULL, init_FF = NULL) {

  warm <- !is.null(init_LL)
  stopifnot(init %in% INITS)

  if (is.infinite(cc)) {
    Ys <- as(Y, "CsparseMatrix")
    if (warm) {
      L0 <- init_LL; F0 <- init_FF
    } else if (init == "rank1") {
      r1 <- fastTopics:::fit_pnmf_rank1(Ys)
      L0 <- cbind(r1$L, matrix(RANK1_PAD, nrow(Ys), K_FIT - 1L))
      F0 <- cbind(r1$F, matrix(RANK1_PAD, ncol(Ys), K_FIT - 1L))
    } else {
      ri <- random_init(nrow(Ys), ncol(Ys), K_FIT)
      L0 <- ri$LL; F0 <- ri$FF
    }
    rownames(L0) <- rownames(Ys); rownames(F0) <- colnames(Ys)
    colnames(L0) <- colnames(F0) <- paste0("k", seq_len(K_FIT))

    fit <- fastTopics::fit_poisson_nmf(
      Ys, fit0 = fastTopics::init_poisson_nmf(Ys, F = F0, L = L0),
      numiter = MAXITER, method = "scd",
      control = list(nc = THREADS, min.delta.loglik = TOL), verbose = "none")

    n_iter <- nrow(fit$progress)
    return(list(LL = fit$L, FF = fit$F, lam = tcrossprod(fit$L, fit$F),
                n_iter = n_iter, converged = n_iter < MAXITER,
                fitter = "fastTopics::fit_poisson_nmf"))
  }

  args <- list(Y = Y, K = K_FIT, cc = cc, s = FALSE, loglik = "exact",
               control = list(maxiter = MAXITER, tol = TOL, max_time = MAX_TIME,
                              verbose = FALSE, threads = THREADS))
  if (warm) { args$init_LL <- init_LL; args$init_FF <- init_FF }
  else      { args$init_method <- init }

  fit <- do.call(log1pNMF::fit_poisson_log1p_nmf, args)
  list(LL = fit$LL, FF = fit$FF, lam = stats::fitted(fit),   # s == 1, so lambda
       n_iter = length(fit$objective_trace) - 1L,
       converged = isTRUE(fit$converged),
       fitter = "log1pNMF::fit_poisson_log1p_nmf")
}

## ---- factor matching + recovery metrics ------------------------------------

#' All permutations of 1:k, as a list. k = 5 => 120, so brute force is fine and
#' avoids a dependency on clue / RcppHungarian.
perms_of <- function(k) {
  if (k == 1L) return(list(1L))
  out <- list()
  for (i in seq_len(k))
    for (rest in perms_of(k - 1L)) {
      r <- rest; r[r >= i] <- r[r >= i] + 1L
      out[[length(out) + 1L]] <- c(i, r)
    }
  out
}

safe_cor <- function(A, B) {
  C <- suppressWarnings(stats::cor(A, B))
  C[!is.finite(C)] <- 0
  C
}

#' Column permutation of `fit` maximizing the total correlation with `truth`.
match_factors <- function(truth, fit) {
  C  <- safe_cor(truth, fit)
  pp <- perms_of(ncol(C))
  sc <- vapply(pp, function(q) sum(C[cbind(seq_len(ncol(C)), q)]), numeric(1))
  pp[[which.max(sc)]]
}

pois_loglik <- function(Y, lam) {
  lam <- pmax(lam, 1e-300)
  sum(Y * log(lam) - lam - lgamma(Y + 1))
}

#' Poisson KL from the true rates to the fitted rates, per entry. This is the
#' parameterization-free measure of how well the mean structure was recovered.
rate_kl <- function(lam_true, lam_fit) {
  lt <- pmax(lam_true, 1e-12); lf <- pmax(lam_fit, 1e-12)
  mean(lt * log(lt / lf) - lt + lf)
}

#' Everything we score a fit on. `perm` is chosen on L (the structure that is
#' shared across every c_true within a seed) and then reused for F.
#'
#' F is scored against the realized truth F_c = alpha*F0 + beta. Correlation is
#' invariant to that affine map, so this is identical to scoring against the
#' sparse F0 -- no separate F0 metric is needed.
score_fit <- function(d, LL_fit, FF_fit, lam_fit) {
  perm <- match_factors(d$truth$L, LL_fit)
  Lc   <- diag(safe_cor(d$truth$L,   LL_fit[, perm, drop = FALSE]))
  Fc   <- diag(safe_cor(d$FF_true,   FF_fit[, perm, drop = FALSE]))

  data.frame(
    L_cor_mean  = mean(Lc),  L_cor_min  = min(Lc),
    F_cor_mean  = mean(Fc),  F_cor_min  = min(Fc),
    rate_kl        = rate_kl(d$lambda, lam_fit),
    rate_cor_log1p = as.numeric(safe_cor(as.vector(log1p(d$lambda)),
                                         as.vector(log1p(lam_fit)))),
    loglik         = pois_loglik(d$Y, lam_fit),
    loglik_oracle  = pois_loglik(d$Y, d$lambda),
    perm           = paste(perm, collapse = "-"))
}

## ---- task id <-> (seed, c_true, c_fit, init) -------------------------------
## 0-based SLURM_ARRAY_TASK_ID over seeds x c_true x c_fit x init, init fastest
## then c_fit. Both initializations of a cell are therefore adjacent task ids,
## so they tend to land on the same node in the same window.
N_TASKS <- N_SEEDS * length(C_GRID) * length(C_GRID) * length(INITS)

task_spec <- function(id) {
  nc <- length(C_GRID); ni <- length(INITS)
  stopifnot(id >= 0, id < N_TASKS)
  init_i <- id %% ni;              rest <- id %/% ni
  cfit_i <- rest %% nc;            rest <- rest %/% nc
  ctru_i <- rest %% nc
  seed_i <- rest %/% nc
  list(seed   = seed_i + 1L,
       c_true = C_GRID[ctru_i + 1L],
       c_fit  = C_GRID[cfit_i + 1L],
       init   = INITS[init_i + 1L])
}

tag_of <- function(spec) {
  fmt <- function(x) if (is.infinite(x)) "Inf" else format(x, scientific = FALSE, trim = TRUE)
  sprintf("cgrid_seed%02d_ctrue%s_cfit%s_%s", spec$seed,
          fmt(spec$c_true), fmt(spec$c_fit), spec$init)
}
