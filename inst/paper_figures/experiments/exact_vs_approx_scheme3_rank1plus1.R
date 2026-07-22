# Exact-vs-approx initialization-timing experiment -- SCHEME 3 of 3:
# rank-1 warm-up + one exact iteration.
#
# Both methods (approximate and exact objective) start from the shared rank-1
# warm-up followed by one exact iteration (the fit_lsa_approximate_supp.R recipe).
# The warm-up and the extra exact step are recorded and prepended to each method's
# from-scratch curve.
#
# One of three independently-launchable scripts; shared config, data loading, and
# helpers live in exact_vs_approx_common.R.
#
# Run from this directory, e.g.:
#   Rscript exact_vs_approx_scheme3_rank1plus1.R > scheme3_rank1plus1.log 2>&1 &

source("exact_vs_approx_common.R")

r1 <- make_rank1_init()

## one additional exact iteration from the padded rank-1 init (tracked)
message("Shared rank-1 + 1 exact iteration")
onestep <- fit_poisson_log1p_nmf(
  Y = Y, cc = cc, loglik = "exact",
  init_LL = r1$rank1_LL, init_FF = r1$rank1_FF, control = ctrl(1)
)

message("Scheme 3: rank-1 + 1 exact iteration initialization")
results <- fit_both(onestep$LL, onestep$FF, "rank1_plus1")
save_scheme("rank1_plus1", results, prefix = list(r1$rank1, onestep))

message("Done (scheme 3: rank1_plus1).")
