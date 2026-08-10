# Compute per-fit summary metrics for the exact-vs-approx simulation and write
# one flat table (RDS + CSV) for downstream plotting. Run on the cluster from
# this directory, after the array is done:
#
#   Rscript analyze_approx_sim.R [out_basename]
#
# Default output: <OUT_DIR>/approx_sim_metrics.{rds,csv} -- one row per fit,
# i.e. per (seed, c_true, c_fit, method); N_TASKS rows (2520 in the extended
# round: 30 seeds x 3 c_true x (9 finite c_fit x 3 methods + Inf, exact
# only)). Small enough to scp home; the full fits (factors, objective
# traces) stay on the server.
#
# Every column the worker stored on `row` is carried through unchanged --
# notably loglik_exact (the EXACT Poisson log-likelihood of the fitted rates,
# whichever objective produced them), loglik_oracle (same under the true
# rates), n_iter, converged, seconds.
#
# ---------------------------------------------------------------------------
# WHAT IS COMPUTED, per fit
#
#   hoyer_L, hoyer_F      mean Hoyer (2004) sparsity over the K fitted
#                         columns: (sqrt(n) - ||x||_1/||x||_2)/(sqrt(n) - 1);
#                         0 = flat, 1 = single nonzero. All-zero columns are
#                         NA and drop from the mean.
#
#   rankcor_L, rankcor_F  mean Spearman correlation over the K(K-1)/2 pairs
#                         of fitted columns -- the redundancy / separation
#                         measure (high = factors carry near-duplicate
#                         rankings).
#
#   L_cor_mean, L_cor_min the recovery measure: fitted columns are matched to
#   F_cor_mean, F_cor_min the stored truth by the column permutation
#                         maximizing total Spearman correlation with L_true
#                         (chosen on L, reused for F -- L is the structure
#                         shared across c_true within a seed), then the
#                         mean/min of the per-column Spearman correlations.
#   perm                  the matching permutation, for spot checks.
#
#   hoyer_L_true, hoyer_F_true, rankcor_L_true, rankcor_F_true
#                         the same sparsity / redundancy statistics for the
#                         stored truth. Constant across (c_fit, method) within
#                         a (seed, c_true); repeated on every row so that
#                         "sparser or denser than the truth" is a subtraction
#                         at plotting time.
#
# Hoyer is scale-invariant per column and Spearman is monotone-invariant, so
# none of this depends on the arbitrary column scalings of L and F.
#
# ---------------------------------------------------------------------------
# SECOND OUTPUT: <OUT_DIR>/approx_sim_bounds.{rds,csv} -- one row per
# (fit, factor k), N_TASKS x K rows, for studying the appendix sparsity
# bounds against the fitted matrices. Per row, with B = L F^t of that fit
# (at c_fit = Inf, B is the fitted rate matrix itself -- the identity-link
# case of the same bounds):
#
#   d_L, d_F         column densities d(x) = ||x||_1 / (sqrt(m) ||x||_2);
#                    H(x) = kappa_m (1 - d(x)), so Hoyer is an affine map of
#                    these (kappa_m = sqrt(m)/(sqrt(m)-1))
#   H_L, H_F         the same columns' Hoyer sparsities, for convenience
#   w                w_k = mean(l_k) mean(f_k) / mean(B); sums to 1 over k
#   pm_L, pm_F       max/mean of the column
#   pm_B             max/mean of B (repeated across the fit's K rows)
#   slack_pm         pm(B) / (pm_L * pm_F * w)      >= 1  (Lemma pm-budget)
#   slack_hyper      d_L * d_F / sqrt(w / pm_B)     >= 1  (Theorem hyperbola)
#
# slack_* near 1 means the corresponding bound is ACTIVE for that column;
# large values mean the fit sits well inside the constraint. The main table
# also gains a pm_B column.

source("approx_sim_common.R")

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1] else file.path(OUT_DIR, "approx_sim_metrics")

## ---- metric helpers (as in analyze_c_grid_sim.R) ---------------------------

hoyer <- function(x) {
  x <- abs(x)
  n <- length(x)
  l2 <- sqrt(sum(x^2))
  if (n < 2L || !is.finite(l2) || l2 == 0) return(NA_real_)
  (sqrt(n) - sum(x) / l2) / (sqrt(n) - 1)
}

mean_hoyer <- function(M) mean(apply(M, 2, hoyer), na.rm = TRUE)

mean_rank_cor <- function(M) {
  if (ncol(M) < 2L) return(NA_real_)
  C <- suppressWarnings(stats::cor(M, method = "spearman"))
  v <- C[upper.tri(C)]
  if (!any(is.finite(v))) return(NA_real_)
  mean(v, na.rm = TRUE)
}

spearman_xcor <- function(A, B) {
  C <- suppressWarnings(stats::cor(A, B, method = "spearman"))
  C[!is.finite(C)] <- 0
  C
}

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

#' Column permutation of `fit` maximizing total Spearman correlation with
#' `truth` (K = 5: 120 permutations, brute force is instant).
match_factors <- function(truth, fit) {
  C  <- spearman_xcor(truth, fit)
  pp <- perms_of(ncol(C))
  sc <- vapply(pp, function(q) sum(C[cbind(seq_len(ncol(C)), q)]), numeric(1))
  pp[[which.max(sc)]]
}

dens  <- function(x) sum(x) / (sqrt(length(x)) * sqrt(sum(x^2)))
pmr   <- function(x) max(x) / mean(x)
kap   <- function(m) sqrt(m) / (sqrt(m) - 1)

#' Per-column bound quantities for one fit (columns in the fit's own order;
#' these are per-column summaries, so no matching to the truth is involved).
bound_rows <- function(LL, FF) {
  B    <- tcrossprod(LL, FF)
  pm_B <- pmr(B)
  mB   <- mean(B)
  data.frame(
    k    = seq_len(ncol(LL)),
    d_L  = apply(LL, 2, dens),  d_F  = apply(FF, 2, dens),
    H_L  = apply(LL, 2, hoyer), H_F  = apply(FF, 2, hoyer),
    w    = colMeans(LL) * colMeans(FF) / mB,
    pm_L = apply(LL, 2, pmr),   pm_F = apply(FF, 2, pmr),
    pm_B = pm_B)
}

## ---- walk the tasks --------------------------------------------------------

specs <- lapply(0:(N_TASKS - 1L), task_spec)
tags  <- vapply(specs, tag_of, character(1))
files <- file.path(OUT_DIR, paste0(tags, ".rds"))

message("reading ", N_TASKS, " fits from ", OUT_DIR)

out         <- vector("list", N_TASKS)
bnd         <- vector("list", N_TASKS)
missing     <- character(0)
truth_cache <- list()          # keyed by "seed|c_true"

for (i in seq_len(N_TASKS)) {
  if (!file.exists(files[i])) { missing <- c(missing, tags[i]); next }
  x <- tryCatch(readRDS(files[i]), error = function(e) NULL)
  if (is.null(x) || is.null(x$row)) { missing <- c(missing, tags[i]); next }

  key <- paste(x$row$seed, format(x$row$c_true), sep = "|")
  if (is.null(truth_cache[[key]]))
    truth_cache[[key]] <- list(
      hoyer_L_true   = mean_hoyer(x$L_true),
      hoyer_F_true   = mean_hoyer(x$F_true),
      rankcor_L_true = mean_rank_cor(x$L_true),
      rankcor_F_true = mean_rank_cor(x$F_true))
  tc <- truth_cache[[key]]

  LL <- as.matrix(x$LL); FF <- as.matrix(x$FF)
  perm <- match_factors(x$L_true, LL)
  Lc <- diag(spearman_xcor(x$L_true, LL[, perm, drop = FALSE]))
  Fc <- diag(spearman_xcor(x$F_true, FF[, perm, drop = FALSE]))

  br <- bound_rows(LL, FF)
  bnd[[i]] <- cbind(x$row[c("task", "seed", "c_true", "cc", "method")], br,
                    row.names = NULL)

  out[[i]] <- cbind(x$row, pm_B = br$pm_B[1], data.frame(
    hoyer_L    = mean_hoyer(LL),
    hoyer_F    = mean_hoyer(FF),
    rankcor_L  = mean_rank_cor(LL),
    rankcor_F  = mean_rank_cor(FF),
    L_cor_mean = mean(Lc), L_cor_min = min(Lc),
    F_cor_mean = mean(Fc), F_cor_min = min(Fc),
    perm       = paste(perm, collapse = "-"),
    hoyer_L_true   = tc$hoyer_L_true,
    hoyer_F_true   = tc$hoyer_F_true,
    rankcor_L_true = tc$rankcor_L_true,
    rankcor_F_true = tc$rankcor_F_true))

  if (i %% 80L == 0L) message("  ", i, " / ", N_TASKS)
}

res <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
rownames(res) <- NULL

if (length(missing))
  message("WARNING: ", length(missing), " fit(s) unreadable or absent: ",
          paste(utils::head(missing, 5), collapse = ", "))

saveRDS(res, paste0(base, ".rds"))
utils::write.csv(res, paste0(base, ".csv"), row.names = FALSE)
message("Wrote ", base, ".rds and ", base, ".csv  (", nrow(res), " fits, ",
        ncol(res), " columns)")

## the per-(fit, k) bound table, with the slack of each appendix inequality
bounds <- do.call(rbind, bnd[!vapply(bnd, is.null, logical(1))])
rownames(bounds) <- NULL
bounds$slack_pm    <- bounds$pm_B / (bounds$pm_L * bounds$pm_F * bounds$w)
bounds$slack_hyper <- bounds$d_L * bounds$d_F / sqrt(bounds$w / bounds$pm_B)
bbase <- file.path(dirname(base), "approx_sim_bounds")
saveRDS(bounds, paste0(bbase, ".rds"))
utils::write.csv(bounds, paste0(bbase, ".csv"), row.names = FALSE)
message("Wrote ", bbase, ".rds and ", bbase, ".csv  (", nrow(bounds),
        " rows = fits x K)")

## ---- a quick look so problems surface before you download ------------------

cat("\nfits per (c_true, c_fit) -- N_SEEDS x 3 methods at finite c_fit,",
    "N_SEEDS at Inf:\n")
print(table(c_true = res$c_true, c_fit = res$cc))

cat("\nEXACT method, mean over seeds (rows = c_true, cols = c_fit):\n")
ex <- res[res$method == "exact", ]
cat("\nhoyer_L:\n");   print(round(tapply(ex$hoyer_L,   list(ex$c_true, ex$cc), mean), 3))
cat("\nhoyer_F:\n");   print(round(tapply(ex$hoyer_F,   list(ex$c_true, ex$cc), mean), 3))
cat("\nrankcor_L:\n"); print(round(tapply(ex$rankcor_L, list(ex$c_true, ex$cc), mean), 3))
cat("\nL_cor_mean:\n"); print(round(tapply(ex$L_cor_mean, list(ex$c_true, ex$cc), mean), 3))

cat("\napproximation gap: mean loglik_exact minus the exact method's, same\n")
cat("(seed, c_true, c_fit); negative = approximation found a worse optimum:\n")
key <- paste(res$seed, res$c_true, res$cc, sep = "|")
ex_ll <- res$loglik_exact[res$method == "exact"][match(key, key[res$method == "exact"])]
gap <- res$loglik_exact - ex_ll
for (m in setdiff(METHODS, "exact")) {
  cat("\n", m, ":\n", sep = "")
  print(signif(tapply(gap[res$method == m],
                      list(res$c_true[res$method == m], res$cc[res$method == m]),
                      mean), 3))
}

cat("\nbound slack sanity: both must be >= 1 everywhere\n")
cat(sprintf("  min slack_pm    = %.4f\n  min slack_hyper = %.4f\n",
            min(bounds$slack_pm), min(bounds$slack_hyper)))
cat("\nmedian slack_hyper, EXACT fits (rows = c_true, cols = c_fit):\n")
be <- bounds[bounds$method == "exact", ]
print(round(tapply(be$slack_hyper, list(be$c_true, be$cc), median), 2))

cat("\nany non-finite metric values?\n")
num <- setdiff(names(res)[vapply(res, is.numeric, logical(1))], c("c_true", "cc"))
bad <- vapply(res[num], function(v) sum(!is.finite(v)), integer(1))
if (any(bad > 0)) print(bad[bad > 0]) else cat("  none\n")
