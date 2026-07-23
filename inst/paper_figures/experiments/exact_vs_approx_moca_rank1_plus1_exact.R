# Exact-vs-approx timing experiment on MOCA.
# Initialization: full rank-1 (rank-1 warm-up + 1 exact iteration).  Objective: exact.
#
# One of two independently-launchable MOCA jobs; shared config/data-loading in
# exact_vs_approx_moca_common.R, shared machinery in exact_vs_approx_helpers.R.
# (Change "rank1_plus1" -> "rank1" below to drop the extra exact init iteration.)
#
# Run from this directory, e.g.:
#   Rscript exact_vs_approx_moca_rank1_plus1_exact.R > moca_rank1_plus1_exact.log 2>&1 &

source("exact_vs_approx_moca_common.R")

run_job("rank1_plus1", "exact")
