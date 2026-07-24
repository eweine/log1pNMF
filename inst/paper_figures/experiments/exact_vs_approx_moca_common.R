# Shared setup for the exact-vs-approx timing experiment on the Mouse
# Organogenesis Cell Atlas (Cao et al. 2019; ~1.33M cells x 26,183 genes).
#
# Same design as exact_vs_approx_common.R (pancreas), but the data arrives as a
# GENES x CELLS dgCMatrix that we transpose to CELLS x GENES on the fly (no
# separate export). Sourced by the MOCA job scripts, which call run_job(...).
#
# MEMORY: the counts matrix is ~10 GB, and the transpose briefly needs a second
# copy (~20 GB peak); the fit then builds triplet summaries + per-row/col
# structures over ~904M nonzeros (tens of GB more). Give the job ~96-128 GB
# (see run_schemes_moca.sbatch). This will NOT run on a 16 GB machine.
#
# EXACT-FIT COST: the exact objective evaluates densely over all 1.33M cells, so
# each exact iteration is far heavier than on the pancreas data. maxiter_exact is
# set small here on purpose; raise it only with a correspondingly long walltime.

library(Matrix)
library(log1pNMF)

## ---- configuration ----------------------------------------------------------
counts_path    <- "/rafalab/eweine/log1pNMF/inst/paper_figures/data/gene_count_cleaned.RDS"  # genes x cells dgCMatrix
K              <- 25        # rank of the factorization (adjust for MOCA as desired)
cc             <- 1         # link-function tuning parameter c
maxiter_approx <- 500       # high safety cap; max_time is the real stop
maxiter_exact  <- 500       # high safety cap; max_time is the real stop
init_maxiter   <- 5         # rank-1 warm-up iterations (package default)
seed           <- 1
n_threads      <- 64
out_dir        <- "/rafalab/eweine/log1p_experiments/"  # MUST exist -- outputs are written here at job end
out_tag        <- sprintf("exact_vs_approx_moca_K%d_c%s", K, format(cc))
tol            <- 1e-8
max_time       <- 10 * 3600 # stop each fit after 10 h of optimization time
                            # (node-speed-invariant; excludes setup + tracking)

## ---- load data: genes x cells  ->  cells x genes ----------------------------
message("Loading MOCA counts from ", counts_path)
counts <- readRDS(counts_path)                 # genes x cells
message(sprintf("  raw: %d genes x %d cells, %s nonzeros",
                nrow(counts), ncol(counts), format(length(counts@x), big.mark = ",")))

## drop all-zero genes (this atlas has one gene with zero counts across all
## cells; such a feature is degenerate for the model). tabulate over the row
## indices avoids a full rowSums pass. All cells have >= 17 nonzeros, so there
## are no all-zero observations to drop.
gkeep <- tabulate(counts@i + 1L, nbins = nrow(counts)) > 0
if (!all(gkeep)) {
  message("  dropping ", sum(!gkeep), " all-zero gene(s)")
  counts <- counts[gkeep, ]
}

message("Transposing to cells x genes (peak memory ~2x here) ...")
Y <- Matrix::t(counts)                          # cells x genes (observations x features)
rm(counts); gc()
n <- nrow(Y)
p <- ncol(Y)
message(sprintf("Data: %d cells x %d genes, %s nonzeros",
                n, p, format(length(Y@x), big.mark = ",")))

config <- list(
  counts_path = counts_path, K = K, cc = cc,
  maxiter_approx = maxiter_approx, maxiter_exact = maxiter_exact,
  init_maxiter = init_maxiter, seed = seed,
  n_threads = n_threads, tol = tol, max_time = max_time, n = n, p = p
)

## shared fitting/init/trace machinery (ctrl, fit_one, run_job, ...)
source("exact_vs_approx_helpers.R")
