# Exact-vs-approx timing experiment on MOCA.
# Initialization: rank-1 warm-up (no extra exact iteration).  Objective: exact.
#
# One of two independently-launchable MOCA jobs; shared config/data-loading in
# exact_vs_approx_moca_common.R, shared machinery in exact_vs_approx_helpers.R.
#
# Submit with run_moca_exact.sbatch.

source("exact_vs_approx_moca_common.R")

run_job("rank1", "exact")
