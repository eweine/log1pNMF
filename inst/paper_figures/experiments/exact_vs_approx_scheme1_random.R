# Exact-vs-approx initialization-timing experiment -- SCHEME 1 of 3: random init.
#
# Both methods (approximate and exact objective) start from a *shared random*
# initialization and are timed while their exact log-likelihood is tracked.
#
# One of three independently-launchable scripts (schemes: random, rank1,
# rank1_plus1) so all three can run in parallel. Shared config, data loading, and
# helpers live in exact_vs_approx_common.R.
#
# Run from this directory, e.g.:
#   Rscript exact_vs_approx_scheme1_random.R > scheme1_random.log 2>&1 &

source("exact_vs_approx_common.R")

## shared random initialization -- drawn exactly as the package's
## init_method = "random" does (runif in [1e-8, 0.05]), generated here so both
## methods receive the same matrices
message("Scheme 1: random initialization")
set.seed(seed)
rand_LL <- matrix(stats::runif(n * K, 1e-8, 0.05), n, K)
rand_FF <- matrix(stats::runif(p * K, 1e-8, 0.05), p, K)

results <- fit_both(rand_LL, rand_FF, "random")
save_scheme("random", results, prefix = list())

message("Done (scheme 1: random).")
