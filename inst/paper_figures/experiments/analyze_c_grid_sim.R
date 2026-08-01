# Compute per-fit summary statistics for every c-grid cell and write them out as
# one flat table (RDS + CSV) for downstream plotting.
#
#   Rscript analyze_c_grid_sim.R [out_basename]
#
# Default output: <OUT_DIR>/c_grid_sim_metrics.{rds,csv}
#
# One row per FIT, i.e. per (seed, c_true, c_fit, init) -- 10 x 8 x 8 x 2 = 1280
# rows. Every column already on the stored `rows` table is carried through, plus
# the metrics below.
#
# ---------------------------------------------------------------------------
# WHAT IS COMPUTED
#
#   loglik            Poisson log-likelihood of Y under the fitted rates.
#                     Already stored by the worker; carried through unchanged.
#                     Directly comparable across c_fit -- same Y, same
#                     dominating measure, constants included.
#
#   rmse              sqrt(mean((lambda_fit - lambda_true)^2)) over all n*p
#                     entries, on the RATE scale.
#   rmse_log1p        the same on the log1p rate scale. Included because the
#                     rates span orders of magnitude, so the raw RMSE is
#                     dominated by the handful of largest entries and can look
#                     flat across c_fit even when the low-rate structure differs
#                     a lot. Use whichever answers the question you are asking.
#
#   hoyer_L, hoyer_F  mean Hoyer sparsity over the K columns of the fitted
#                     loadings / factors. Hoyer (2004) for a length-n vector:
#                         (sqrt(n) - ||x||_1/||x||_2) / (sqrt(n) - 1)
#                     0 when all entries are equal, 1 when one entry carries
#                     everything. An all-zero column is undefined and returns NA
#                     (dropped from the mean).
#
#   rankcor_F         mean Spearman correlation between the K fitted factor
#                     columns, averaged over the K*(K-1)/2 distinct pairs (3
#                     pairs at K = 3). High values mean the recovered factors
#                     are redundant / not well separated.
#   rankcor_L         the same for the loadings. Not requested, but it is the
#                     exact analogue and free to compute; ignore the column if
#                     it is not useful.
#
#   *_true            the same sparsity / rank-correlation statistics for the
#                     TRUE L and F_c. Constant across c_fit and init within a
#                     (seed, c_true), but repeated on every row so that
#                     "did the fit come out sparser or denser than the truth"
#                     is a single subtraction at plotting time.
#
# ---------------------------------------------------------------------------
# NOTES
#
# lambda_true is regenerated from (seed, c_true) rather than stored -- the
# generator is deterministic, so this reproduces the exact rates the counts were
# drawn from. The calibration solve costs ~1 s, so it is done once per
# (seed, c_true) and reused across that pair's 8 c_fit x 2 init fits.
#
# lambda_fit is reconstructed from the saved LL / FF, matching what the fitter
# would have returned: cc * expm1(LL F' / max(1, cc)) for finite c_fit, and
# LL F' for c_fit = Inf (the linear limit, where fastTopics fits the rate
# directly). No size factors anywhere, so lambda is the mean.

suppressMessages(library(Matrix))
source("c_grid_sim_common.R")

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1] else file.path(OUT_DIR, "c_grid_sim_metrics")

## ---- metric helpers --------------------------------------------------------

#' Hoyer sparsity of a vector: 0 = flat, 1 = a single nonzero entry.
hoyer <- function(x) {
  x <- abs(x)
  n <- length(x)
  l2 <- sqrt(sum(x^2))
  if (n < 2L || !is.finite(l2) || l2 == 0) return(NA_real_)
  (sqrt(n) - sum(x) / l2) / (sqrt(n) - 1)
}

mean_hoyer <- function(M) mean(apply(M, 2, hoyer), na.rm = TRUE)

#' Mean Spearman correlation over the distinct column pairs.
mean_rank_cor <- function(M) {
  if (ncol(M) < 2L) return(NA_real_)
  C <- suppressWarnings(stats::cor(M, method = "spearman"))
  v <- C[upper.tri(C)]
  if (!any(is.finite(v))) return(NA_real_)
  mean(v, na.rm = TRUE)
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

#' Fitted rates from the saved factors. Mirrors fitted.log1p_nmf_fit (with
#' s == 1) for finite c, and the plain product in the c -> Inf limit.
lambda_from_fit <- function(LL, FF, cc) {
  th <- tcrossprod(LL, FF)
  if (is.infinite(cc)) return(th)
  cc * expm1(th / max(1, cc))
}

## ---- walk the grid ---------------------------------------------------------

specs <- lapply(seq_len(N_TASKS) - 1L, task_spec)
tags  <- vapply(specs, tag_of, character(1))
files <- file.path(OUT_DIR, paste0(tags, ".rds"))

message("reading ", N_TASKS, " cells from ", OUT_DIR)

out     <- vector("list", N_TASKS)
missing <- character(0)
truth_cache <- list()      # keyed by "seed|c_true"; holds lambda + truth stats

for (i in seq_len(N_TASKS)) {
  sp <- specs[[i]]
  if (!file.exists(files[i])) { missing <- c(missing, tags[i]); next }
  x <- tryCatch(readRDS(files[i]), error = function(e) NULL)
  if (is.null(x) || is.null(x$rows) || is.null(x$fits)) {
    missing <- c(missing, tags[i]); next
  }

  ## truth for this (seed, c_true), built once and reused across c_fit / init
  key <- paste(sp$seed, sp$c_true, sep = "|")
  if (is.null(truth_cache[[key]])) {
    d <- sim_dataset(sp$seed, sp$c_true)
    truth_cache[[key]] <- list(
      lambda       = d$lambda,
      lambda_log1p = log1p(d$lambda),
      hoyer_L_true   = mean_hoyer(d$truth$L),
      hoyer_F_true   = mean_hoyer(d$FF_true),
      hoyer_F0_true  = mean_hoyer(d$truth$F0),
      rankcor_L_true = mean_rank_cor(d$truth$L),
      rankcor_F_true = mean_rank_cor(d$FF_true))
  }
  tc <- truth_cache[[key]]

  rows <- x$rows
  add  <- do.call(rbind, lapply(seq_len(nrow(rows)), function(r) {
    ini <- rows$init[r]
    f   <- x$fits[[ini]]
    if (is.null(f)) return(NULL)
    lam <- lambda_from_fit(f$LL, f$FF, sp$c_fit)
    data.frame(
      rmse        = rmse(lam, tc$lambda),
      rmse_log1p  = rmse(log1p(pmax(lam, 0)), tc$lambda_log1p),
      hoyer_L     = mean_hoyer(f$LL),
      hoyer_F     = mean_hoyer(f$FF),
      rankcor_L   = mean_rank_cor(f$LL),
      rankcor_F   = mean_rank_cor(f$FF),
      hoyer_L_true   = tc$hoyer_L_true,
      hoyer_F_true   = tc$hoyer_F_true,
      hoyer_F0_true  = tc$hoyer_F0_true,
      rankcor_L_true = tc$rankcor_L_true,
      rankcor_F_true = tc$rankcor_F_true)
  }))

  out[[i]] <- cbind(rows, add)
  if (i %% 80L == 0L) message("  ", i, " / ", N_TASKS, " cells")
}

res <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
rownames(res) <- NULL

if (length(missing))
  message("WARNING: ", length(missing), " cell(s) unreadable or absent, e.g. ",
          paste(utils::head(missing, 3), collapse = ", "))

saveRDS(res, paste0(base, ".rds"))
utils::write.csv(res, paste0(base, ".csv"), row.names = FALSE)
message("Wrote ", base, ".rds and ", base, ".csv  (", nrow(res), " fits, ",
        ncol(res), " columns)")

## ---- a quick look so problems surface before you download ------------------
cat("\ncolumns:\n"); print(names(res))
cat("\nfits per (c_true, c_fit) -- should be 20 (10 seeds x 2 inits):\n")
print(table(c_true = res$c_true, c_fit = res$c_fit))
cat("\nmean rmse_log1p (rows = c_true, cols = c_fit), best init per cell:\n")
b <- res[order(paste(res$seed, res$c_true, res$c_fit), -res$loglik), ]
b <- b[!duplicated(paste(b$seed, b$c_true, b$c_fit)), ]
print(signif(tapply(b$rmse_log1p, list(b$c_true, b$c_fit), mean), 3))
## c_true / c_fit legitimately contain Inf, so exclude them or this cries wolf
cat("\nsample of the computed metrics:\n")
print(head(res[order(res$c_true, res$c_fit, res$init),
               c("c_true","c_fit","init","loglik","rmse","rmse_log1p",
                 "hoyer_L","hoyer_F","rankcor_F","hoyer_F_true","rankcor_F_true")],
           8), row.names = FALSE, digits = 4)

cat("\nany non-finite metric values?\n")
num  <- setdiff(names(res)[vapply(res, is.numeric, logical(1))], c("c_true", "c_fit"))
bad  <- vapply(res[num], function(v) sum(!is.finite(v)), integer(1))
if (any(bad > 0)) print(bad[bad > 0]) else cat("  none\n")
