# Exact-vs-approx initialization-timing experiment.
# Scheme: rank1 initialization.  Objective: approx.
#
# One of six independently-launchable jobs (3 schemes x 2 objectives) so every
# fit runs in its own job and is checkpointed independently. Shared config, data
# loading, and helpers live in exact_vs_approx_common.R.
#
# Run from this directory, e.g.:
#   Rscript exact_vs_approx_rank1_approx.R > rank1_approx.log 2>&1 &

source("exact_vs_approx_common.R")

run_job("rank1", "approx")
