# Exact-vs-approx initialization-timing experiment -- SCHEME 2 of 3: rank-1 init.
#
# Both methods (approximate and exact objective) start from a *shared rank-1
# warm-up* (the package's rank1 procedure, padded to rank K). The warm-up
# trajectory is recorded and prepended to each method's from-scratch curve.
#
# One of three independently-launchable scripts; shared config, data loading, and
# helpers live in exact_vs_approx_common.R.
#
# Run from this directory, e.g.:
#   Rscript exact_vs_approx_scheme2_rank1.R > scheme2_rank1.log 2>&1 &

source("exact_vs_approx_common.R")

r1 <- make_rank1_init()

message("Scheme 2: rank-1 initialization")
results <- fit_both(r1$rank1_LL, r1$rank1_FF, "rank1")
save_scheme("rank1", results, prefix = list(r1$rank1))

message("Done (scheme 2: rank1).")
